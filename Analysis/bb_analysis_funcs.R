
# Returns bootstratpped correlations
boot_cor = function(x, y, ci = .95, nsamples = 10000) {
  if(length(x) != length(y)) {
    stop("Vectors must be of equal length")
  }
  N = length(x)
  samps = mclapply(1:nsamples,
                   function(i) {
                     idxs = sample(1:N, N, replace=T)
                     this.x = x[idxs]
                     this.y = y[idxs]
                     return(cor(this.x, this.y))
                   })
  samps = unlist(samps)
  min_quant = (1-ci)/2
  max_quant = ci + (1-ci)/2
  return(quantile(samps,c(min_quant,max_quant)))
}


# A helper function that turns off legends (for nicer plotting through Illustrator)
no_legend = theme(legend.position='none',strip.background = element_blank(),strip.text = element_blank())

# Adds a 1:1 line
abl = geom_abline(intercept = 0, slope = 1, linetype = 'dashed')

# Dodge
pd = position_dodge(.1)

# Saves figures easily without a legend
save_figure = function(figure, filename, width = 4, height = 4, remove_legend=T) {
  ofl = paste(FIGURE_DIR, filename, sep='')
  if(remove_legend) {figure = figure + no_legend}
  ggsave(ofl, figure, width = width, height = height, units = 'in')
}

# Plot accuracy graphs
plot_raw_acc = function(dat, xaxis = NULL, color = NULL, use_emp = T, xlab = NULL) {
  
  rdat = ungroup(dat)
  rdat$N = with(rdat, NLeft + NBal + NRight)
  accnm = ifelse(use_emp, 'EmpAcc','ModAcc')
  rdat = select_(rdat, xaxis, color, accnm, 'N')
  names(rdat) = c('XAxis','Color','Acc','N')
  gdat = rdat %>% group_by(XAxis, Color) %>%
    summarize(Acc = mean(Acc), N = sum(N)) %>%
    mutate(SE = sqrt(Acc*(1-Acc)/N))
  
  plt = ggplot(data=gdat, mapping=aes(x=XAxis, y=Acc, color=Color, group=Color)) +
    geom_line(position=pd) + geom_pointrange(aes(ymin = Acc-SE, ymax = Acc+SE),position=pd) + 
    scale_y_continuous(labels = scales::percent,limits = c(0,1)) +
    xlab(xlab) + ylab('Accuracy') + theme_bw() +
    scale_color_discrete(drop = F)
  
  return(plt)
}

plot_acc_cor = function(dat, color = NULL, shape = NULL, xlab = 'Model Accuracy', ylab = 'Empirical Accuracy',
                        alpha=1.) {
  return(ggplot(data=dat, mapping = aes_string(x="ModAcc", y="EmpAcc", color=color, shape=shape)) +
           abl + geom_point(alpha=alpha) + 
           scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
           scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
           xlab(xlab) + ylab(ylab) + 
           theme_bw() + scale_color_discrete(drop=F))
}

# Plot choice graphs
plot_choices = function(dat, xpanel = NULL, ypanel = NULL, ptsize = 1, ylab = NULL) {
  
  # Relevel
  if(length(factor(dat$Type)) == 6) {
    dat$Type = factor(dat$Type, c("Bal", "Weight", "Dist", "CB", "CW", "CD"))
  }
  
  
  # Calculate the model predictions by group
  mdat = dat %>% select_(xpanel, ypanel, 'ModLeft', 'ModBal', 'ModRight')
  names(mdat) = c('Xpan','Ypan','L','B','R')
  mdat = mdat %>% group_by(Xpan,Ypan) %>% 
    summarize(L = mean(L), B = mean(B), R = mean(R)) %>%
    gather(Choice, Model, -Xpan, -Ypan)
  
  # Cacluate empirical predictions by group
  edat = dat %>% mutate(N = NLeft + NBal + NRight) %>% 
    select_(xpanel, ypanel, 'EmpLeft', 'EmpBal', 'EmpRight', 'N')
  names(edat) = c('Xpan','Ypan','L','B','R', 'N')
  edat = edat %>% group_by(Xpan,Ypan) %>%
    summarize(L = mean(L), B = mean(B), R = mean(R), N = sum(N)) %>%
    gather(Choice, Empirical, -Xpan, -Ypan, -N)
  
  gdat = merge(edat, mdat)
  gdat$EmpSE = with(gdat, sqrt(Empirical*(1-Empirical)/N))
  gdat$Choice = factor(as.character(gdat$Choice), levels = c('L','B','R'))
  
  plt = ggplot(data=gdat, mapping = aes(x=Choice, y = Empirical)) +
    geom_bar(stat='identity') + 
    geom_linerange(aes(ymin = Empirical - EmpSE, ymax = Empirical + EmpSE)) +
    geom_point(aes(y = Model), color='red', size = ptsize, alpha = .7) +
    facet_grid(Ypan ~ Xpan) + theme_bw() + ylab(ylab) +
    scale_y_continuous(breaks=c(0,.5,1), limits=c(-.02,1))
  
  return(plt)
  
}

plot_choices_ind = function(dat, wid, xpanel = NULL, ypanel = NULL, ptsize = 1, ylab = NULL) {
  dat %>% filter(WID == wid) %>% mutate(NLeft = EmpLeft, NBal = EmpBal, NRight = EmpRight) %>%
    plot_choices(xpanel, ypanel, ptsize, ylab)
}

# Return a table with LLH comparisons for each of the different model combinations

# NOTE: MAY NEED TO STRIP OUT PARAMETERS THAT ARE MOSTLY UNUSED
get_n_params = function(modtype = 'shapes', subtype = 'smp') {
  shps = c('m' = 8, 's' = 7, 'p' = 8, 'sm' = 9, 'ms' = 9, 
           'mp' = 10, 'pm' = 10, 'ps' = 10, 'sp' = 10)
  matps = shps + 3
  balps = c('m' = 11, 's' = 7, 'p' = 12, 'sm' = 12, 'ms' = 12, 
            'mp' = 14, 'pm' = 14, 'ps' = 14, 'sp' = 14)
  if(modtype == 'shapes') {
    if (nchar(subtype) == 3) {
      return(11)
    }
    
    return(shps[subtype])
  }
  
  if(modtype == 'materials') {
    if(nchar(subtype) == 3) {
      return(13)
    }
    return(matps[subtype])
  }
  
  if(modtype == 'balance') {
    if(nchar(subtype) == 3) {
      return(15)
    }
    return(balps[subtype])
  }
}

modcomp_llh_table = function(folder, prefix, nsamps = 5, modtype = 'shapes', do_ci = T) {
  
  if(modtype == 'shapes') {
    loadfn = loadShapesAgg
  } else if(modtype == 'materials') {
    loadfn = loadMatAgg
  } else if(modtype == 'balance') {
    loadfn = loadBalAgg
  } else if(modtype == 'combined') {
    loadfn = loadCombAgg
  }
  
  modtypes = c('smp','spm','sm','sp','s','msp','mps','ms','mp','m','psm','pms','ps','pm','p')
  tab = NULL
  for(mod in modtypes) {
    pflnm = paste(folder,'/',prefix,'_',mod,'_',nsamps,'_params.csv',sep='')
    llh = -read.csv(pflnm)$LLH
    thisdat = loadfn(paste(folder,'/',prefix,'_',mod,'_',nsamps,'.csv',sep=''))
    bic = -2*llh + get_n_params(modtype,mod)*log(nrow(thisdat))
    
    if(do_ci) {
      bootcors = (with(thisdat,boot_cor(EmpAcc, ModAcc, ci=.95, nsamples=10000)))
      tab = rbind(tab, data.frame(Model = mod, LLH = llh, BIC = bic,
                                  r = with(thisdat, cor(EmpAcc,ModAcc)), r_0.025 = bootcors[1], r_0.975 = bootcors[2]))
    } else {
      tab = rbind(tab, data.frame(Model = mod, LLH = llh, BIC = bic, r = with(thisdat, cor(EmpAcc,ModAcc))))
    }
    
  }
  rownames(tab) = NULL
  return(tab)
}

modcomp_across_table = function(folder_prefix, file_prefix, nsamps = c(1,2,5,10), modtype = 'shapes') {
  
  df = modcomp_llh_table(paste(folder_prefix,nsamps[1],sep='_'), file_prefix, nsamps[1], modtype, F)
  names(df)[2:3] = c(paste('Samp',nsamps[1],sep='_'), paste('BIC_Samp',nsamps[1],sep='_'))
  bic_df = df[c(1,3)]
  df = df[1:2]
  
  for (s in nsamps[2:length(nsamps)]) {
    ndf = modcomp_llh_table(paste(folder_prefix,s,sep='_'), file_prefix, s, modtype, F)
    names(ndf)[2:3] = c(paste('Samp',s,sep='_'), paste('BIC_Samp',s,sep='_'))
    bic_df = merge(bic_df, ndf[c(1,3)])
    df = merge(df, ndf[1:2])
  }
  bestllh = with(df, sapply(1:nrow(df), function(i) {return(max(df[i,2:ncol(df)]))}))
  bestidx = with(df, sapply(1:nrow(df), function(i) {return(which(df[i,2:ncol(df)] == bestllh[i]))}))
  df$BestSamps = names(df)[bestidx + 1]
  df$BestLLH = bestllh
  df$BestBIC = sapply(1:nrow(df), function(i) {return(bic_df[i,bestidx[i]+1])})
  
  return(df)
}

### NOTE::: NEED TO FIX UP WITH COMPLETELY CONTROLLED TRIALS
modcomp_cv_multisamples = function(folder_prefix, nsims = c(1,2,3,4,5,10)) {
  
  # Get the baseline cv-llhs
  bfl = read.csv(paste(folder_prefix,baseline[2],'sim/cv_llh_grid.csv',sep=''))
  
  
  for(ns in nsims) {
    
    fl = read.csv(paste(folder_prefix,ns,'sim/cv_llh_grid.csv',sep=''))
    new_llh = fl %>% select(-Order) %>% sapply(mean)
    
    if(isFirst) {
      best_llh = new_llh
      best_sims = rep(ns, times=length(best_llh))
      isFirst = F
    } else {
      bs_idx = which(new_llh > best_llh)
      best_sims[bs_idx] = ns
      best_llh = mapply(max, best_llh, new_llh)
    }
    
  }
  
  r = data.frame(Model = names(best_llh), AvgCVLLH = best_llh, NSims = best_sims) %>% arrange(desc(AvgCVLLH))
  return(r)
}

modcomp_cv_bysims = function(folder_prefix, nsims = c(1,2,3,4,5,10)) {
  
  r = NULL
  for (ns in nsims) {
    fl = read.csv(paste(folder_prefix,ns,'sim/cv_llh_grid.csv',sep=''))
    new_llh = fl %>% select(-Order) %>% sapply(mean)
    r = rbind(r, c(Sims=ns,new_llh))
  }
  return(data.frame(r))
}

modcomp_cv_full = function(folder_prefix, nsims = c(1,2,3,4,5,10)) {
  
  r = NULL
  for(ns in nsims) {
    fl = read.csv(paste(folder_prefix,ns,'sim/cv_llh_grid.csv',sep=''))
    r = rbind(r, fl %>% mutate(Sims = ns))
  }
  return(r)
}

apply_modcomp_base = function(fulldat, basemod = 'msp', basesims = 1) {
  
  baseline = subset(fulldat, Sims == basesims)[[basemod]]
  
  trans_dat = cbind(fulldat %>% select(Order, Sims),
                    sapply(fulldat %>% select(-Order, -Sims), function(x) {x-baseline})) %>%
    gather(Model, LLH, -Order, -Sims)
  
  avg_llhs = trans_dat %>% group_by(Sims, Model) %>% summarize(AvgLLH = mean(LLH)) %>% spread(Model, AvgLLH) %>% ungroup
  sims = avg_llhs$Sims
  nms = names(avg_llhs %>% select(-Sims))
  best_sims = sims[avg_llhs %>% select(-Sims) %>% sapply(function(x) {which(x == max(x))})]
  
  best_dat = subset(trans_dat, paste(Model,Sims) %in% paste(nms, best_sims))
  best_dat$Sims = factor(best_dat$Sims)
  return(best_dat)
}

plot_cv_modtype = function(llhtab, sort_llhtab = NULL) {
  if(is.null(sort_llhtab)) {sort_llhtab = llhtab}
  # Set up the plot data frame
  avgllhs = sort_llhtab %>% select(-Order) %>% sapply(mean) %>% sort(decreasing = T)
  d = llhtab %>% select(-Order) %>% gather(key=Model, value=LLH)
  d$Model = factor(d$Model, levels = names(avgllhs))
  
  pl = ggplot(d, aes(x = Model, y = LLH)) + geom_violin() + theme_bw() + ylab('Cross-validated LLH')
  return(pl)
}

plot_oversims_cv_modtype = function(folder_prefix, nsims = c(1,2,3,4,5,10), modtype = 'smp') {
  
  p = NULL
  for (ns in nsims) {
    
    fl = read.csv(paste(folder_prefix,ns,'sim/cv_llh_grid.csv',sep=''))
    p = rbind(p, data.frame(Sims = ns, Order = fl[["Order"]], CV_LLH = fl[[modtype]]))
    
  }
  p$Sims = factor(p$Sims)
  p$Order = factor(p$Order)

  ggplot(p, aes(x = Order, y = CV_LLH, group = Sims, color=Sims)) + geom_line() + theme_bw() + ylim(c(-7500, -5500))
  
}

plot_cv_addind = function(llhtab) {
  base = llhtab$Baseline
  llhdiff = llhtab %>% select(-Order, -Baseline) %>% lapply(function(x) {x-base}) %>% data.frame
  
  d = llhdiff %>% gather(Parameter, LLH)
  
  pl = ggplot(d, aes(x=Parameter, y = LLH)) + geom_hline(yintercept = 0) + geom_violin() + 
    theme_bw() + ylab('Cross-validated LLH')
  return(pl)
}

boot_acc_noiseceil_cor = function(dat, nsamps = 100) {
  
  uwids = unique(dat$WID)
  N = length(uwids)
  
  spl_cor = function() {
    ud = dat %>% mutate(InS1 = WID %in% sample(uwids, floor(N/2), replace=F)) %>%
      group_by(Trial, InS1) %>% summarize(Acc = mean(EmpAcc)) %>%
      spread(InS1, Acc)
    return(cor(ud['FALSE'], ud['TRUE']))
  }
  
  cors = replicate(nsamps, spl_cor())
  return(mean(cors))
  
}

boot_acc_model_asceil = function(dat_ind, dat_agg, nsamps = 100) {
  
  uwids = unique(dat_ind$WID)
  N = length(uwids)
  
  spl_cor = function() {
    ud = dat_ind %>% filter(WID %in% sample(uwids, floor(N/2), replace=F)) %>%
      group_by(Trial) %>% summarize(Acc = mean(EmpAcc)) %>%
      merge(dat_agg %>% select(Trial, ModAcc))
    return(with(ud, cor(Acc, ModAcc)))
  }
  
  cors = replicate(nsamps, spl_cor())
  return(mean(cors))
  
}

boot_compare_vs_model = function(dat_ind, dat_agg, nsamps = 100) {
  uwids = unique(dat_ind$WID)
  N = length(uwids)
  
  spl_cor = function() {
    ud = dat_ind %>% mutate(InS1 = ifelse(WID %in% sample(uwids, floor(N/2), replace=F), 'S1','S2')) %>%
      group_by(Trial,InS1) %>% summarize(Acc = mean(EmpAcc)) %>%
      spread(InS1, Acc) %>%
      merge(dat_agg %>% select(Trial, ModAcc))
    return(with(ud, c(cor(S1, S2), cor(S1, ModAcc))))
  }
  
  cors = replicate(nsamps, spl_cor())
  rets = apply(cors, 1, mean)
  names(rets) = c('SplitHalf','VsMod')
  return(rets)
}

boot_mse_vs_model = function(dat_ind, dat_agg, nsamps = 100) {
  uwids = unique(dat_ind$WID)
  N = length(uwids)
  
  spl_cor = function() {
    ud = dat_ind %>% mutate(InS1 = ifelse(WID %in% sample(uwids, floor(N/2), replace=F), 'S1','S2')) %>%
      group_by(Trial,InS1) %>% summarize(Acc = mean(EmpAcc)) %>%
      spread(InS1, Acc) %>%
      merge(dat_agg %>% select(Trial, ModAcc))
    return(with(ud, c(sum((S1-S2)^2), sum((S1-ModAcc)^2))))
  }
  
  cors = replicate(nsamps, spl_cor())
  rets = apply(cors, 1, mean)
  names(rets) = c('SplitHalf','VsMod')
  return(rets)
}


boot_compare_vs_rules_and_model = function(dat_ind, dat_agg, rules_agg, nsamps = 100) {
  uwids = unique(dat_ind$WID)
  N = length(uwids)
  
  spl_cor = function() {
    ud = dat_ind %>% mutate(InS1 = ifelse(WID %in% sample(uwids, floor(N/2), replace=F), 'S1','S2')) %>%
      group_by(Trial,InS1) %>% summarize(Acc = mean(EmpAcc)) %>%
      spread(InS1, Acc) %>%
      merge(dat_agg %>% select(Trial, ModAcc)) %>%
      merge(rules_agg %>% mutate(RuleAcc = ModAcc) %>% select(Trial, RuleAcc))
    return(with(ud, c(cor(S1, S2), cor(S1, ModAcc), cor(S1, RuleAcc))))
  }
  
  cors = replicate(nsamps, spl_cor())
  rets = apply(cors, 1, mean)
  names(rets) = c('SplitHalf','VsMod', 'VsRules')
  return(rets)
}

boot_mse_rules_and_model = function(dat_ind, dat_agg, rules_agg, nsamps = 100) {
  uwids = unique(dat_ind$WID)
  N = length(uwids)
  
  spl_cor = function() {
    ud = dat_ind %>% mutate(InS1 = ifelse(WID %in% sample(uwids, floor(N/2), replace=F), 'S1','S2')) %>%
      group_by(Trial,InS1) %>% summarize(Acc = mean(EmpAcc)) %>%
      spread(InS1, Acc) %>%
      merge(dat_agg %>% select(Trial, ModAcc)) %>%
      merge(rules_agg %>% mutate(RuleAcc = ModAcc) %>% select(Trial, RuleAcc))
    return(with(ud, c(sum((S1-S2)^2), sum((S1-ModAcc)^2), sum((S1-RuleAcc)^2))))
  }
  
  cors = replicate(nsamps, spl_cor())
  rets = apply(cors, 1, mean)
  names(rets) = c('SplitHalf','VsMod', 'VsRules')
  return(rets)
}

make_boot_noiseceil_df = function(dat, group_var, nsamps = 100) {
  group = dat[[group_var]]
  vs = unique(group)
  bcs = sapply(vs, function(v) {boot_acc_noiseceil_cor(dat[group==v,], nsamps)})
  r = data.frame(a = vs, Noise_Ceiling = bcs)
  names(r)[1] = group_var
  return(r)
}

plot_vs_rules = function(basic_dat, rules = c('R1','R2','R3','R4'), incl_rulemod = F, strict_off = F) {
  
  types = c('Bal','Weight','Dist','CW','CD','CB')
  
  # Make the settled rules
  rule_df = data.frame(Type = factor(types, levels=types),
                       R1 = c(1,1,0,1,0,0),
                       R2 = c(1,1,1,1,0,0),
                       R3 = c(1,1,1,1/3,1/3,1/3),
                       R4 = c(1,1,1,1,1,1)) %>% gather(Rule,RuleAcc,-Type)
  
  rule_left = data.frame(Type = types,
                         R1 = c(0,1,0,1,0,1),
                         R2 = c(0,1,1,1,0,1),
                         R3 = c(0,1,1,1/3,1/3,1/3),
                         R4 = c(0,1,1,1,1,0)) %>% gather(Rule,RuleLeft,-Type)
  
  rule_bal = data.frame(Type = types,
                        R1 = c(1,0,1,0,0,0),
                        R2 = c(1,0,0,0,0,0),
                        R3 = c(1,0,0,1/3,1/3,1/3),
                        R4 = c(1,0,0,0,0,1)) %>% gather(Rule,RuleBal,-Type)
  
  rule_right = data.frame(Type = types,
                          R1 = c(0,0,0,0,1,0),
                          R2 = c(0,0,0,0,1,0),
                          R3 = c(0,0,0,1/3,1/3,1/3),
                          R4 = c(0,0,0,0,0,0)) %>% gather(Rule,RuleRight,-Type)
  
  rules_agg = rule_left %>% merge(rule_bal) %>% merge(rule_right)
  
  broad_dat = basic_dat %>% group_by(Type, Rule) %>% summarize(EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight),
                                                               ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight), 
                                                               RMLeft = mean(RM_Left), RMBal = mean(RM_Bal), RMRight = mean(RM_Right), N = length(WID)) %>% 
    mutate(EmpLSD = sqrt(EmpLeft*(1-EmpLeft)/N), EmpBSD = sqrt(EmpBal*(1-EmpBal)/N), EmpRSD = sqrt(EmpRight*(1-EmpRight))/N)  %>% merge(rules_agg)
  
  pl_dat = broad_dat %>% mutate(Choice = 'L', Pct = EmpLeft, SD = EmpLSD, RuleCh = RuleLeft, ModCh = ModLeft, RMCh = RMLeft) %>% select(Type, Rule, Choice, Pct, SD, RuleCh, ModCh, RMCh) %>%
    rbind(broad_dat %>% mutate(Choice = 'B', Pct = EmpBal, SD = EmpBSD, RuleCh = RuleBal, ModCh = ModBal, RMCh = RMBal) %>% select(Type, Rule, Choice, Pct, SD, RuleCh, ModCh, RMCh)) %>% 
    rbind(broad_dat %>% mutate(Choice = 'R', Pct = EmpRight, SD = EmpRSD, RuleCh = RuleRight, ModCh = ModRight, RMCh = RMRight) %>% select(Type, Rule, Choice, Pct, SD, RuleCh, ModCh, RMCh))
  
  
  pl_dat$Choice = factor(pl_dat$Choice, levels=c('L','B','R'))
  
  r = ggplot(pl_dat, aes(x = Choice, y = Pct, ymin=Pct-2*SD, ymax = Pct+2*SD)) + geom_bar(stat='identity') + geom_linerange() + 
    geom_point(aes(y=ModCh), color='red', size=2, alpha = .8) + 
    facet_grid(Type~ Rule) + theme_bw()
  if(!strict_off) {
    r = r + geom_point(aes(y=RuleCh), color='blue', size=2, alpha = .8)
  }
  if(incl_rulemod) {
    r = r + geom_point(aes(y=RMCh), color='green', size=2, alpha = .8)
  }
  return(r)
  
}

plot_vs_pquartile = function(basic_dat, masspdat) {
  udat = basic_dat %>% select(-Rule) %>% 
    merge(masspdat %>% mutate(Rule = Mass_Prob_Quantile) %>% select(WID, Rule))
  
  broad_dat = udat %>% group_by(Type, Rule) %>% summarize(EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight),
                                                               ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight), 
                                                               RMLeft = mean(RM_Left), RMBal = mean(RM_Bal), RMRight = mean(RM_Right), N = length(WID)) %>% 
    mutate(EmpLSD = sqrt(EmpLeft*(1-EmpLeft)/N), EmpBSD = sqrt(EmpBal*(1-EmpBal)/N), EmpRSD = sqrt(EmpRight*(1-EmpRight))/N)
  
  pl_dat = broad_dat %>% mutate(Choice = 'L', Pct = EmpLeft, SD = EmpLSD, ModCh = ModLeft, RMCh = RMLeft) %>% select(Type, Rule, Choice, Pct, SD, ModCh, RMCh) %>%
    rbind(broad_dat %>% mutate(Choice = 'B', Pct = EmpBal, SD = EmpBSD, ModCh = ModBal, RMCh = RMBal) %>% select(Type, Rule, Choice, Pct, SD, ModCh, RMCh)) %>% 
    rbind(broad_dat %>% mutate(Choice = 'R', Pct = EmpRight, SD = EmpRSD, ModCh = ModRight, RMCh = RMRight) %>% select(Type, Rule, Choice, Pct, SD, ModCh, RMCh))
  
  
  pl_dat$Choice = factor(pl_dat$Choice, levels=c('L','B','R'))
  
  r = ggplot(pl_dat, aes(x = Choice, y = Pct, ymin=Pct-2*SD, ymax = Pct+2*SD)) + geom_bar(stat='identity') + geom_linerange() + 
    geom_point(aes(y=ModCh), color='red', size=2, alpha = .8) + geom_point(aes(y=RMCh), color='green', size=2, alpha = .8) + 
    facet_grid(Type~ Rule) + theme_bw()
  return(r)
}

plot_full_v_rulesonly = function(fulldat, rulesdat, xpanel, ypanel) {
  mdat = fulldat %>% merge(rulesdat %>% mutate(RuleLeft = ModLeft, RuleBal = ModBal,
                                               RuleRight = ModRight) %>%
                             select(WID, Trial, Rule, RuleLeft, RuleBal, RuleRight))
  
  broad_dat = mdat %>% group_by_(xpanel, ypanel) %>% summarize(EmpLeft = mean(EmpLeft), EmpBal = mean(EmpBal), EmpRight = mean(EmpRight),
                                                          ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight), 
                                                          RuleLeft = mean(RuleLeft), RuleBal = mean(RuleBal), RuleRight = mean(RuleRight), N = length(WID)) %>% 
    mutate(EmpLSD = sqrt(EmpLeft*(1-EmpLeft)/N), EmpBSD = sqrt(EmpBal*(1-EmpBal)/N), EmpRSD = sqrt(EmpRight*(1-EmpRight))/N)
 
  pl_dat = broad_dat %>% mutate(Choice = 'L', Pct = EmpLeft, SD = EmpLSD, ModCh = ModLeft, RuleCh = RuleLeft) %>% select_(xpanel, ypanel, "Choice", "Pct", "SD", "ModCh", "RuleCh") %>%
    rbind(broad_dat %>% mutate(Choice = 'B', Pct = EmpBal, SD = EmpBSD, ModCh = ModBal, RuleCh = RuleBal) %>% select_(xpanel, ypanel, "Choice", "Pct", "SD", "ModCh", "RuleCh")) %>%
    rbind(broad_dat %>% mutate(Choice = 'R', Pct = EmpRight, SD = EmpRSD, ModCh = ModRight, RuleCh = RuleRight) %>% select_(xpanel, ypanel, "Choice", "Pct", "SD", "ModCh", "RuleCh"))
  
  names(pl_dat)[1:2] = c('xpan','ypan')
  pl_dat$Choice = factor(pl_dat$Choice, levels=c('L','B','R'))
  
  r = ggplot(pl_dat, aes(x = Choice, y = Pct, ymin=Pct-2*SD, ymax = Pct+2*SD)) + geom_bar(stat='identity') + geom_linerange() + 
    geom_point(aes(y=ModCh), color='red', size=2, alpha = .8) + geom_point(aes(y=RuleCh), color='green', size=2, alpha = .8) + 
    facet_grid(xpan ~ ypan) + theme_bw()
  return(r)
  
}


# Function for plotting individual parameter cross-validations

plot_indparam = function(dat, lower = .1, upper = .9) {
  
  all_params = c('distance_uncert', 'tower_weight_uncert', 'beam_weight', 'pivot_range', 'pivot_range_2.5',
                 'pivot_range_5', 'pivot_range_10', 'pivot_range_20', 'weight_jnd', 'distance_jnd', 'p_weight',
                 'phys_uncert', 'lapse_rate')
  
  gdat = dat %>% gather(param, dllh, -Order) %>% group_by(param) %>% 
    summarize(avg = mean(dllh), bottom = quantile(dllh, lower), top = quantile(dllh, upper))
  
  these_params = all_params[all_params %in% gdat$param]
  
  gdat$param = factor(gdat$param, levels=these_params)
  
  ggplot(gdat, aes(x=param, y=avg, ymin = bottom, ymax = top)) + geom_hline(yintercept=0) + geom_linerange() + 
    geom_point() + theme_bw() + xlab('Parameters') + ylab(expression(paste(Delta, "LLH")))
}
