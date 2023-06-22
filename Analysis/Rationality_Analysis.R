library(parallel)
library(lme4)
library(ggplot2)
library(tidyr)
library(dplyr)
library(xtable) # Not directly used but helpful for latex translating
library(data.table)
library(RColorBrewer)
library(jsonlite)


# Many of the functions for loading trials are hidden in the `loadData.R' file -- 
# these transform the raw data from the model into formats used for this analysis
source('load_bb_data.R')

rat_dat = read.csv("../Modeling/output/rationality/mod_predictions.csv") %>% 
  filter(Experiment != "combined")

params = read_json('../Modeling/output/comb_strats/base_strats_all_params.json')

FIT_WEIGHT_HEURISTIC = with(params$strategies, smp / (smp + sp))

# Figure out the correct answer for these figures
correct_falls = loadShapesAgg('../Modeling/output/comb_strats/base_strats_all_sh.csv') %>% select(Trial, Falls, Type) %>% mutate(Experiment='shapes') %>%
  rbind(loadMatAgg('../Modeling/output/comb_strats/base_strats_all_mat.csv') %>% select(Trial, Falls, Type) %>% mutate(Experiment='materials')) %>% 
  rbind(loadBalAgg('../Modeling/output/comb_strats/base_strats_all_bal.csv') %>% select(Trial, Falls, Type) %>% mutate(Experiment='pivot')) %>%
  rbind(loadCombAgg('../Modeling/output/extensions_from_comb/base_strats_combined.csv') %>% select(Trial, Falls, Type) %>% mutate(Experiment='combined'))

rat_dat = merge(rat_dat, correct_falls)


# Calculate the VoC for each trial
value_of_computation = function(dat, sym_cost, mass_cost, dist_cost, phys_cost) {
  dat$VoC = with(dat, ifelse(Falls=='B', ModBal, ifelse(Falls=='L', ModLeft, ModRight)) -
                   (UsesSym * sym_cost + UsesMass * mass_cost + UsesDist * dist_cost + UsesPhys * phys_cost))
  return(dat)
}

# Calculate the rate of return for each trial
rate_of_return = function(dat, sym_time, mass_time, dist_time, phys_time, inter_trial_time) {
  dat = dat %>%
    mutate(E_Points = ifelse(Falls=='B', ModBal, ifelse(Falls=='L', ModLeft, ModRight)),
           CalcTime = UsesSym*sym_time + UsesMass*mass_time + UsesDist*dist_time + UsesPhys*phys_time + inter_trial_time,
           RoR = E_Points / CalcTime)
  return(dat)
}


rat_cost = value_of_computation(rat_dat %>% filter(!grepl("d", OrderType)),
                                .05, .2, .2, .3)
rat_cost %>% group_by(OrderType) %>% summarize(Utility = mean(VoC))

######################
# Overall accuracy

# Physics gets it right ~52% of the time
rat_dat %>% filter(OrderType == 'p') %>% 
  mutate(Acc = ifelse(Falls=='B', ModBal, ifelse(Falls=='L', ModLeft, ModRight))) %>%
  with(summary(Acc))

# Symmetry tops out at 16%
# Weight tops out at 1/3rd (and fails rather than pass through in 1/6th)

#######################
# Make plots

get_full_voc_smp = function(dat, sym_cost, mass_cost, phys_cost) {
  calc = value_of_computation(dat %>% filter(!grepl("d", OrderType)), sym_cost, mass_cost, 0, phys_cost) %>% group_by(OrderType) %>% 
           summarize(Utility=mean(VoC))
  return(with(calc, as.list(structure(Utility, names=as.character(OrderType)))))
}

basic_voc_grid = data.frame(sym_cost = rep(seq(.05, .5, by=.05), each=100),
                            mass_cost = rep(seq(.05, .5, by=.05), times=10, each=10),
                            phys_cost = rep(seq(.05, .5, by=.05), times=100))

basic_voc_grid = as.data.table(basic_voc_grid)[, get_full_voc_smp(rat_dat, 
                                                                  sym_cost, mass_cost, phys_cost), 
                                               by= .(sym_cost, mass_cost, phys_cost)]

basic_voc_grid$BestOrder = colnames(basic_voc_grid)[4:ncol(basic_voc_grid)][max.col(basic_voc_grid[,4:ncol(basic_voc_grid)], ties.method='first')]

approx_eq = function(a,b,tol=1e-8) {
  return(abs(a-b)<tol)
}
make_voc_grid_by_phys = function(p_cost, dat = basic_voc_grid) {
  ggplot(dat %>% filter(approx_eq(phys_cost, p_cost)), aes(x=mass_cost, y=sym_cost, fill=BestOrder)) +
    geom_tile()
  
}


#+ ########################
# Analyze by sets of trials

get_voc_trial_smp = function(dat, sym_cost, mass_cost, phys_cost) {
  calc = value_of_computation(dat %>% filter(!grepl("d", OrderType)), sym_cost, mass_cost, 0, phys_cost) %>%
    select(Experiment, Trial, OrderType, VoC) %>% spread(OrderType, VoC)
  calc$BestOrder = colnames(calc)[3:ncol(calc)][max.col(calc[,3:ncol(calc)], ties.method='first')]
  return(calc)
}

get_distrib_voc_smp = function(dat, sym_cost, mass_cost, phys_cost) {
  return(with(get_voc_trial_smp(dat, sym_cost, mass_cost, phys_cost),
              table(Experiment, BestOrder)))
}

get_resampled_winners_smp = function(dat, sym_cost, mass_cost, phys_cost, N=100) {
  calc = value_of_computation(dat %>% filter(!grepl("d", OrderType)), sym_cost, mass_cost, 0, phys_cost) %>%
    select(Experiment, Trial, OrderType, VoC) %>% spread(OrderType, VoC)
  
  do_resamp = function(exp) {
    if (exp != 'all') {
      sdat = calc %>% filter(Experiment == exp)
    } else {
      sdat = calc
    }
    n_rows = nrow(sdat)
    rs_dat = sdat[sample.int(n_rows, replace=T),]
    tots = rs_dat %>% select(-Experiment, -Trial) %>% summarise_each(funs(sum))
    return(colnames(tots)[max.col(tots, ties.method = 'first')])
  }
  
  ret_dat = data.frame(Experiment = rep(c('shapes','materials','pivot','combined','all'), each=N),
                       BestSamp = c(replicate(N, do_resamp('shapes')),
                                    replicate(N, do_resamp('materials')),
                                    replicate(N, do_resamp('pivot')),
                                    replicate(N, do_resamp('combined')),
                                    replicate(N, do_resamp('all'))))
  return(with(ret_dat, table(Experiment, BestSamp)))
}

get_resampled_winners_all = function(dat, sym_cost, mass_cost, phys_cost, N=100) {
  calc = value_of_computation(dat %>% filter(!grepl("d", OrderType)), sym_cost, mass_cost, 0, phys_cost) %>%
    select(Experiment, Trial, OrderType, VoC) 
  ots = unique(as.character(calc$OrderType))
  calc = calc %>% spread(OrderType, VoC)
  
  do_resamp = function() {
    sdat = calc
    sdat$guessing = 1/3
    n_rows = nrow(sdat)
    rs_dat = sdat[sample.int(n_rows, replace=T),]
    tots = rs_dat %>% select(-Experiment, -Trial) %>% summarise_each(funs(sum))
    return(colnames(tots)[max.col(tots, ties.method = 'first')])
  }
  counts = structure(rep(0,length(ots)+1), names=c(ots, 'guessing'))
  for(i in 1:N) {
    rs = do_resamp()
    counts[rs] = counts[rs] + 1
  }
  return(as.list(counts))
}


# NOTE: THIS TAKES 10+ MINUTES
#grid_choice_rs_winner = data.frame(sym_cost = rep(seq(.05, .5, by=.05), each=100),
#                            mass_cost = rep(seq(.05, .5, by=.05), times=10, each=10),
#                            phys_cost = rep(seq(.05, .5, by=.05), times=100))

#grid_choice_rs_winner = data.frame(sym_cost = rep(seq(.05, .5, by=.05), each=100),
#                            mass_cost = rep(seq(.05, .5, by=.05), times=10, each=10),
#                            phys_cost = rep(seq(.05, .5, by=.05), times=100))

# grid_choice_rs_winner = data.frame(sym_cost = rep(seq(.02, .2, by=.02), each=100),
#                                    mass_cost = rep(seq(.02, .2, by=.02), times=10, each=10),
#                                    phys_cost = rep(seq(.02, .2, by=.02), times=100))
# 
# 
# grid_choice_rs_winner = as.data.table(grid_choice_rs_winner)[, get_resampled_winners_all(rat_dat, 
#                                                                                 sym_cost, mass_cost, phys_cost), 
#                                                              by= .(sym_cost, mass_cost, phys_cost)]


#+ ########################################
# Look at strategies based on costs

find_reasonable_strategies = function(calculated, thresh = 10) {
  return(calculated %>%
    mutate(strat_msp_smp = (msp + smp) >= thresh,
           strat_sp = sp >= thresh,
           strat_p_only = p >= thresh,
           strat_ms_sm_s_m = (m+s+ms+sm) >= thresh,
           strat_guess = guessing >= thresh)
  )
}

make_costs_from_phys_cost = function(phys_cost, reasonable_thresh = 10) {
  params = data.frame(sym_cost_raw = sequence(1:10)/10,
                      mass_cost_raw = rep(1:10, times=1:10)/10) %>%
    mutate(sym_cost = sym_cost_raw*phys_cost,
           mass_cost = mass_cost_raw*phys_cost,
           phys_cost = phys_cost)
  
  rs_wins = as.data.table(params)[, get_resampled_winners_all(rat_dat,
                                                                             sym_cost, mass_cost, phys_cost),
                                                 by= .(sym_cost, mass_cost, phys_cost, sym_cost_raw, mass_cost_raw)]
  calc_wins = find_reasonable_strategies(rs_wins, reasonable_thresh)
  return(calc_wins)
}

plot_strat_from_phys_costs = function(phys_cost, reasonable_thresh = 10) {
  strats = make_costs_from_phys_cost(phys_cost, reasonable_thresh)
  plt_strat = strats %>% select(sym_cost, mass_cost, strat_msp_smp, strat_sp, strat_p_only, strat_ms_sm_s_m) %>%
    gather(key=Strategy, value=IsGood, -sym_cost, -mass_cost)
  unit_diff = phys_cost / 20
  plt_strat = plt_strat %>%
    mutate(sym_plot = sym_cost - ifelse(Strategy %in% c('strat_p_only', 'strat_ms_sm_s_m'), unit_diff, 0),
           mass_plot = mass_cost - ifelse(Strategy %in% c('strat_msp_smp', 'strat_ms_sm_s_m'), unit_diff, 0),
           strat_fill = as.factor(ifelse(IsGood, as.character(Strategy), "NA")))
  plt = ggplot(plt_strat, aes(x = mass_plot, y=sym_plot, fill=strat_fill)) + geom_tile() +
    xlab("Weight Rule Cost") + ylab("Symmetry Rule Cost") +
    scale_fill_manual(values=c('strat_msp_smp' = "red", 'strat_sp' = "blue", 'strat_p_only' = "darkgreen", 'strat_ms_sm_s_m' = "yellow", "NA" = "white"),
                      labels=c('strat_msp_smp' = "SWP/WSP", 'strat_sp' = "SP", 'strat_p_only' = "Phys Only", 'strat_ms_sm_s_m' = "No Phys", "NA" = "NA")) +
    theme_bw()
  return(list(data=strats, plot=plt))
}

multi_plot_phys_strats = function(phys_costs = seq(.03, .27, by=.03), reasonable_thresh = 10, data=rat_dat) {
  nphyscosts = length(phys_costs)
  all_params = data.frame(phys_cost = rep(phys_costs, each=55)) %>%
    mutate(sym_cost = phys_cost*rep(sequence(1:10)/10, times=nphyscosts),
           mass_cost = phys_cost*rep(rep(1:10, times=1:10)/10, times=nphyscosts))
  rs_wins = as.data.table(all_params)[, get_resampled_winners_all(data,
                                                              sym_cost, mass_cost, phys_cost),
                                  by= .(sym_cost, mass_cost, phys_cost)]
  calc_wins = find_reasonable_strategies(rs_wins, reasonable_thresh)
  plt_strat = calc_wins %>% select(sym_cost, mass_cost, phys_cost, strat_msp_smp, strat_sp, strat_p_only, strat_ms_sm_s_m, strat_guess) %>%
    gather(key=Strategy, value=IsGood, -sym_cost, -mass_cost, -phys_cost) %>%
    mutate(unit_diff = phys_cost /20,
           sym_plot = sym_cost - ifelse(Strategy %in% c('strat_p_only', 'strat_ms_sm_s_m', 'strat_guess'), unit_diff, 0),
           mass_plot = mass_cost - ifelse(Strategy %in% c('strat_msp_smp', 'strat_ms_sm_s_m', 'strat_guess'), unit_diff, 0),
           strat_fill = as.factor(ifelse(IsGood, as.character(Strategy), "NA")),
           box_size = ifelse(Strategy == 'strat_guess', unit_diff*2, unit_diff))
  plt_boxes = calc_wins %>% select(sym_cost, mass_cost, phys_cost) %>%
    mutate(unit_diff = phys_cost / 20,
           bottom = sym_cost - unit_diff,
           left = mass_cost - unit_diff)
  plt = ggplot() + 
    geom_rect(data=plt_strat %>% filter(Strategy=='strat_guess'), aes(xmin = mass_plot, xmax=mass_plot+box_size, ymin=sym_plot, 
                                  ymax=sym_plot+box_size, fill=strat_fill)) +
    geom_rect(data=plt_strat %>% filter(Strategy != 'strat_guess', strat_fill != "NA"), aes(xmin = mass_plot, xmax=mass_plot+box_size, ymin=sym_plot, 
                                  ymax=sym_plot+box_size, fill=strat_fill)) +
    geom_rect(data=plt_boxes, aes(xmin=left, ymin=bottom, xmax=left+unit_diff*2, ymax=bottom+unit_diff*2), 
               alpha=0, fill='white', color='black') +
    xlab("Weight Rule Cost") + ylab("Symmetry Rule Cost") +
    scale_fill_manual(values=c('strat_msp_smp' = "red", 'strat_sp' = "blue", 'strat_p_only' = "darkgreen", 'strat_ms_sm_s_m' = "yellow", "NA" = "white", 'strat_guess' = 'black'),
                      labels=c('strat_msp_smp' = "SWP/WSP", 'strat_sp' = "SP", 'strat_p_only' = "Phys Only", 'strat_ms_sm_s_m' = "No Phys", "NA" = "NA", 'strat_guess' = "Guess")) +
    theme_bw() + facet_wrap(~phys_cost, scales='free')
  return(list(data=calc_wins, plot=plt))
}

bratstrat_fl = "../Modeling/output/rationality/basic_strats.RData"
if(file.exists(bratstrat_fl)) {
  load(bratstrat_fl)
} else {
  strats_for_all = multi_plot_phys_strats()
  strats_for_simple = multi_plot_phys_strats(data = rat_dat %>% filter(Type %in% c("Bal","Dist","Weight")))
  save(strats_for_all, strats_for_simple, file=bratstrat_fl)
}


#+ ##### Do analysis with distance

get_resampled_winners_with_dist = function(dat, sym_cost, mass_cost, dist_cost, phys_cost, N=100) {
  calc = value_of_computation(dat, sym_cost, mass_cost, dist_cost, phys_cost) %>%
    select(Experiment, Trial, OrderType, VoC) 
  ots = unique(as.character(calc$OrderType))
  calc = calc %>% spread(OrderType, VoC)
  do_resamp = function() {
    sdat = calc
    n_rows = nrow(sdat)
    rs_dat = sdat[sample.int(n_rows, replace=T),]
    tots = rs_dat %>% select(-Experiment, -Trial) %>% summarise_each(funs(sum))
    return(colnames(tots)[max.col(tots, ties.method = 'first')])
  }
  counts = structure(rep(0,length(ots)), names=ots)
  for(i in 1:N) {
    rs = do_resamp()
    counts[rs] = counts[rs] + 1
  }
  return(as.list(counts))
}
  
make_costs_for_dist = function(phys_cost, dat = rat_dat, reasonable_thresh = 10) {
  params = data.frame(sym_cost_raw = sequence(1:10)/10,
                      md_cost_raw = rep(1:10, times=1:10)/10) %>%
    mutate(sym_cost = sym_cost_raw*phys_cost,
           mass_dist_cost = md_cost_raw*phys_cost,
           phys_cost = phys_cost)
  
  rs_wins = as.data.table(params)[, get_resampled_winners_with_dist(dat,
                                                              sym_cost, mass_dist_cost, mass_dist_cost, phys_cost),
                                  by= .(sym_cost, mass_dist_cost, phys_cost)]
  #calc_wins = find_reasonable_strategies(rs_wins, reasonable_thresh)
  return(calc_wins)
}

#+ Do analysis by rate ####################

# 
get_pointrate_trial_smp = function(dat, sym_time, mass_time, phys_time, inter_trial_time = 1.) {
  calc = rate_of_return(dat %>% filter(!grepl("d", OrderType)), sym_time, mass_time, 0, phys_time, inter_trial_time) %>%
    select(Experiment, Trial, OrderType, E_Points, CalcTime, RoR)
  return(calc)
}

get_resampled_pointrate_smp = function(dat, sym_time, mass_time, phys_time, inter_trial_time = 1., N=100) {
  calc = rate_of_return(dat %>% filter(!grepl("d", OrderType)), sym_time, mass_time, 0, phys_time, inter_trial_time) %>%
    select(Experiment, Trial, OrderType, E_Points, CalcTime, RoR)
  
  strats = unique(calc$OrderType)
  trialnames = unique(calc$Trial)
  
  split_dat = lapply(strats, function(s) {calc %>% filter(OrderType == s)})
  names(split_dat) = strats
  
  do_one_resamp = function() {
    usedidxs = sample.int(length(trialnames), replace = T)
    rates = sapply(strats, function(s) {
      sdat = split_dat[[as.character(s)]][usedidxs,]
      return(with(sdat, sum(E_Points) / sum(CalcTime)))
    })
    rates = c(rates, (1/3) / inter_trial_time)
    winidx = which(rates == max(rates))[1]
    if (winidx == length(rates)) {
      return('guessing')
    } else {
      return(strats[winidx])
    }
  }
  
  ret = structure(rep(0,length(strats)+1), names=c(as.character(strats), 'guessing'))
  for(i in 1:N) {
    rs = as.character(do_one_resamp())
    ret[rs] = ret[rs] + 1
  }
  return(as.list(ret))
}

make_ror_from_phys_time = function(phys_time, reasonable_thresh = 10) {
  params = data.frame(sym_time = sequence(1:10)/10 * phys_time,
                      mass_time = rep(1:10, times=1:10)/10 * phys_time,
                      phys_time = phys_time)
  
  rs_wins = as.data.table(params)[, get_resampled_pointrate_smp(rat_dat,
                                                              sym_time, mass_time, phys_time),
                                  by= .(sym_time, mass_time, phys_time)]
  rscalc_wins = find_reasonable_strategies(rs_wins, reasonable_thresh)
  return(calc_wins)
}

multi_plot_phys_strats_ror = function(phys_times = seq(.2, 1., by=.1), reasonable_thresh = 10, data=rat_dat) {
  nphyscosts = length(phys_times)
  data = data %>% filter(!grepl("d", OrderType))
  all_params = data.frame(phys_time = rep(phys_times, each=55)) %>%
    mutate(sym_time = phys_time*rep(sequence(1:10)/10, times=nphyscosts),
           mass_time = phys_time*rep(rep(1:10, times=1:10)/10, times=nphyscosts))
  rs_wins = as.data.table(all_params)[, get_resampled_pointrate_smp(data,
                                                                  sym_time, mass_time, phys_time),
                                      by= .(sym_time, mass_time, phys_time)]
  calc_wins = find_reasonable_strategies(rs_wins, reasonable_thresh)
  plt_strat = calc_wins %>% select(sym_time, mass_time, phys_time, strat_msp_smp, strat_sp, strat_p_only, strat_ms_sm_s_m, strat_guess) %>%
    gather(key=Strategy, value=IsGood, -sym_time, -mass_time, -phys_time) %>%
    mutate(unit_diff = phys_time /20,
           sym_plot = sym_time - ifelse(Strategy %in% c('strat_p_only', 'strat_ms_sm_s_m', 'strat_guess'), unit_diff, 0),
           mass_plot = mass_time - ifelse(Strategy %in% c('strat_msp_smp', 'strat_ms_sm_s_m', 'strat_guess'), unit_diff, 0),
           strat_fill = as.factor(ifelse(IsGood, as.character(Strategy), "NA")),
           box_size = ifelse(Strategy == 'strat_guess', unit_diff*2, unit_diff))
  plt_boxes = calc_wins %>% select(sym_time, mass_time, phys_time) %>%
    mutate(unit_diff = phys_time / 20,
           bottom = sym_time - unit_diff,
           left = mass_time - unit_diff)
  plt = ggplot() + 
    geom_rect(data=plt_strat %>% filter(Strategy=='strat_guess'), aes(xmin = mass_plot, xmax=mass_plot+box_size, ymin=sym_plot, 
                                                                      ymax=sym_plot+box_size, fill=strat_fill)) +
    geom_rect(data=plt_strat %>% filter(Strategy != 'strat_guess', strat_fill != "NA"), aes(xmin = mass_plot, xmax=mass_plot+box_size, ymin=sym_plot, 
                                                                                            ymax=sym_plot+box_size, fill=strat_fill)) +
    geom_rect(data=plt_boxes, aes(xmin=left, ymin=bottom, xmax=left+unit_diff*2, ymax=bottom+unit_diff*2), 
              alpha=0, fill='white', color='black') +
    xlab("Weight Rule Time") + ylab("Symmetry Rule Time") +
    scale_fill_manual(values=c('strat_msp_smp' = "red", 'strat_sp' = "blue", 'strat_p_only' = "darkgreen", 'strat_ms_sm_s_m' = "yellow", "NA" = "white", 'strat_guess' = 'black'),
                      labels=c('strat_msp_smp' = "MSP/SMP", 'strat_sp' = "SP", 'strat_p_only' = "Phys Only", 'strat_ms_sm_s_m' = "No Phys", "NA" = "NA", 'strat_guess' = "Guess")) +
    theme_bw() + facet_wrap(~phys_time, scales='free')
  return(list(data=calc_wins, plot=plt))
}

dratstrat_fl = "../Modeling/output/rationality/rate_strats.RData"
if(file.exists(dratstrat_fl)) {
  load(dratstrat_fl)
} else {
  strats_rate_for_all = multi_plot_phys_strats_ror()
  strats_rate_for_simple = multi_plot_phys_strats_ror(data = rat_dat %>% filter(Type %in% c("Bal","Dist","Weight")))
  save(strats_rate_for_all, strats_rate_for_simple, file=(dratstrat_fl))
}



# ##########################
# Assuming correct strats, what's the optimal mix of sp / smp?

calc_optim_sp_smp_mix = function(sym_cost, mass_cost, phys_cost, exp_filter=NULL) {
  if (!is.null(exp_filter)) {
    d = rat_dat %>% filter(Experiment %in% exp_filter)
  } else {
    d = rat_dat
  }
  d = value_of_computation(d %>% filter(OrderType %in% c('sp', 'smp')), sym_cost, mass_cost, 99999999, phys_cost) %>%
    select(Experiment, Trial, OrderType, VoC) %>% spread(OrderType, VoC)
  return(with(d, mean(smp > sp)))
}

stratmix_fl = "../Modeling/output/rationality/mixes.RData"
if(file.exists(stratmix_fl)) {
  load(stratmix_fl)
} else {
  sp_smp_mix = data.frame(phys_cost = rep(seq(.03, .27, by=.03), each=1275)) %>%
    mutate(sym_cost = phys_cost*rep(sequence(1:50)/50, times=9),
           mass_cost = phys_cost*rep(rep(1:50, times=1:50)/50, times=9))
  sp_smp_mix$SMP_Mix = with(sp_smp_mix, mapply(calc_optim_sp_smp_mix,sym_cost, mass_cost, phys_cost))
  sp_smp_mix$Shapes_Mix =with(sp_smp_mix, mapply(calc_optim_sp_smp_mix,sym_cost, mass_cost, phys_cost, 'shapes'))
  sp_smp_mix$Materials_Mix =with(sp_smp_mix, mapply(calc_optim_sp_smp_mix,sym_cost, mass_cost, phys_cost, 'materials'))
  sp_smp_mix$Pivot_Mix =with(sp_smp_mix, mapply(calc_optim_sp_smp_mix,sym_cost, mass_cost, phys_cost, 'pivot'))
  save(sp_smp_mix, file=(stratmix_fl))
}


multiplot_mix = function(column_name = "SMP_Mix") {
  low_cl = colorRampPalette(c('blue', 'red'))
  mid_cl = "#FF0000"
  high_cl = colorRampPalette(c('red', 'yellow'))
  #myscale = scale_fill_gradient2(midpoint=.37, low='blue', mid='red', high='yellow', limits=c(.2, .7))
  # Scale so red is .35 - .45
  myscale = scale_fill_gradientn(colors=c(colorRampPalette(c('blue', 'purple'))(11),  rep("#FF0000",10), colorRampPalette(c('orange', 'yellow'))(24)),  limits=c(.25, .7))
  plt = ggplot(sp_smp_mix, aes_string(x="mass_cost", y="sym_cost", fill=column_name, z=column_name)) +
    geom_raster() +
    facet_wrap(~phys_cost, scales='free') +
    theme_bw() +
    xlab("Weight Rule Cost") +
    ylab("Symmetry Rule Cost") +
    myscale
  return(plt)
}


PLT_DIR = "../Figures/SubFigs/"
ggsave(paste(PLT_DIR, "rat_optimstrats.pdf", sep=''), strats_for_all$plot, units='in', width=10, height=8)
ggsave(paste(PLT_DIR, "rat_optimmix.pdf", sep=''), multiplot_mix(), units='in', width=10, height=8)
