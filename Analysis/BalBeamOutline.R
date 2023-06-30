
#' ---
#' title: Balance Beams -- Outline 
#' author: Kevin A Smith
#' output:
#'    html_document:
#'      toc: true
#'      toc_depth: 1
#'      toc_float: false
#'      theme: default
#'      highlight: tango
#' ---

#+ General settings, echo = FALSE, results = 'hide', fig.width = 4, fig.height = 4 ------------------------------------------------------------------------------

knitr::opts_chunk$set(warning=F, message=F, cache = F, echo=F)
options(digits = 3)
kable = knitr::kable
export = F
FIGURE_DIR = "../Figures/SubFigs/"

#+ Initialization ----------------------------------------------------------

library(parallel)
library(lme4)
library(car)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(xtable) # Not directly used but helpful for latex translating
library(jsonlite)

# Many of the functions for loading trials are hidden in the `loadData.R' file -- 
# these transform the raw data from the model into formats used for this analysis
source('load_bb_data.R')

# Load in the various models used for the Shapes experiment
# * _agg: fit with aggregate model
# * _ind: fit with individual strategy weights
# * _ind_full: with individual strategy weights, but detailed by individual
shapes_agg = loadShapesAgg('../Modeling/anonymized_output/comb_strats/base_strats_all_sh.csv')
shapes_ind = loadShapesInd('../Modeling/anonymized_output/comb_strats/base_strats_individual_sh.csv')
shapes_ind_full = loadShapesInd('../Modeling/anonymized_output/comb_strats/base_strats_individual_sh.csv',T)

materials_agg = loadMatAgg('../Modeling/anonymized_output/comb_strats/base_strats_all_mat.csv')
materials_ind = loadMatInd('../Modeling/anonymized_output/comb_strats/base_strats_individual_mat.csv')
materials_ind_full = loadMatInd('../Modeling/anonymized_output/comb_strats/base_strats_individual_mat.csv', T)

balance_agg = loadBalAgg('../Modeling/anonymized_output/comb_strats/base_strats_all_bal.csv')
balance_ind = loadBalInd('../Modeling/anonymized_output/comb_strats/base_strats_individual_bal.csv')
balance_ind_full = loadBalInd('../Modeling/anonymized_output/comb_strats/base_strats_individual_bal.csv', T)

# Fit using rule methodology - mixture of rules in aggregate
# _ind_: fit allowing a single rule per person
rules_shapes = loadShapesAgg('../Modeling/anonymized_output/comb_strats/rules_all_sh.csv')
rules_materials = loadMatAgg('../Modeling/anonymized_output/comb_strats/rules_all_mat.csv')
rules_balance = loadBalAgg('../Modeling/anonymized_output/comb_strats/rules_all_bal.csv')
rules_ind_shapes = loadShapesInd('../Modeling/anonymized_output/comb_strats/rules_individual_sh.csv',T)
rules_ind_mat = loadMatInd('../Modeling/anonymized_output/comb_strats/rules_individual_mat.csv',T)
rules_ind_bal = loadBalInd('../Modeling/anonymized_output/comb_strats/rules_individual_bal.csv',T)

# Load in best fit parameters by individual
complete_params = read.csv('../Modeling/anonymized_output/ind_strat_choice.csv')
rules_ind_shapes = rules_ind_shapes %>% merge(complete_params %>% select(WID, Rule=rules_Rule))
rules_ind_mat = rules_ind_mat %>% merge(complete_params %>% select(WID, Rule=rules_Rule))
rules_ind_bal = rules_ind_bal %>% merge(complete_params %>% select(WID, Rule=rules_Rule))

# For the geometric experiment
# _geom: assumes no weight information used
geomat_agg = loadMatAgg('../Modeling/anonymized_output/extensions_from_comb/base_strats_geomat.csv', use_geomat = T)
geomat_geom = loadMatAgg('../Modeling/anonymized_output/geomat/base_strats_all_geomat.csv', use_geomat = T)
geomat_raw = read.csv('../Modeling/anonymized_output/raw_data/BB_GeoMatData.csv')

# For the combination experiment
comb_agg = loadCombAgg('../Modeling/anonymized_output/extensions_from_comb/base_strats_combined.csv')
rules_comb = loadCombAgg('../Modeling/anonymized_output/extensions_from_comb/rules_combined.csv')

# Load in crossvalidation data (running many split-halfs)
crossval_dat = read.csv("../Modeling/anonymized_output/strat_crossval/all_crossval.csv")
crossval_ind_vs_rules = read.csv("../Modeling/anonymized_output/strat_crossval/ind_strat_cv_comp.csv") %>% 
  mutate(dLLH = ModLLH_CV - RuleLLH_CV, dPhysLLH = ModLLH_CV - PhysOnlyLLH_CV) %>% 
  merge(complete_params %>% select(WID, Experiment, Rule=rules_Rule))
# crossval_ind_vs_rules = read.csv("../Modeling/output/comb_strats/mod_vs_rules_llh.csv") %>% 
#   spread(Strategy, CVLLH) %>% 
#   mutate(dLLH = base_strats - rules) %>% 
#   merge(complete_params %>% select(WID, Experiment, Rule=rules_Rule))

# Load in data for the Ferretti torque-difference replication
ferretti_agg = loadFerrettiAgg('../Modeling/anonymized_output/extensions_from_comb/base_strats_ferretti.csv')
ferretti_raw = read.csv('../Modeling/anonymized_output/raw_data/BB_FerrettiData.csv')
rules_ferretti = loadFerrettiAgg('../Modeling/anonymized_output/extensions_from_comb/rules_ferretti.csv')

remove_rev = function(tnm) {
  spl = strsplit(tnm, '_')[[1]]
  spl = spl[1:(length(spl)-1)]
  return(paste(spl,collapse='_'))
}
comb_agg_full = read.csv('../Modeling/anonymized_output/raw_data/BB_CombData.csv') %>%
  mutate(EmpLeft = 1*(NormResp=='L'), EmpBal = 1*(NormResp=='B'), EmpRight = 1*(NormResp=='R'),
         EmpAcc = 1*(as.character(NormResp)==as.character(CorrectResp)))
caf_n = table(comb_agg_full$WID)
good_caf = names(caf_n)[caf_n == max(caf_n)]
comb_agg_full$Trial = sapply(as.character(comb_agg_full$TrialName), remove_rev)
comb_agg_full = comb_agg_full %>% subset(WID %in% good_caf) %>%
  select(WID, Trial, EmpLeft, EmpBal, EmpRight, EmpAcc) %>%
  merge(comb_agg %>% select(-EmpLeft, -EmpBal, -EmpRight, -EmpAcc))


# Load in data for shifting strategy experiment & filter:
#  a) People who don't have all trials recorded
#  b) People who responded too quickly (< 1s on over half of trials, indicative of "clicking through")
shift_raw_inp = read.csv("../Modeling/anonymized_output/raw_data/BB_BeneData.csv") %>% 
  mutate(Prioritizes = ifelse(TrialType == 'CD', ifelse(NormResp=='L', 'Dist', ifelse(NormResp=='R', 'Weight', 'Bal')),
                              ifelse(NormResp=='L', 'Weight', ifelse(NormResp=='R', 'Dist', 'Bal'))))

filter_shift_sum = shift_raw_inp %>% group_by(WID) %>% 
  summarize(N = length(WID), MeanRT = mean(Time), MedRT = median(Time))

nrow(filter_shift_sum %>% filter(N < 180))
filter_shift_sum = filter_shift_sum %>% filter(N == 180) # Remove incomplete records first

filter_shift_sum %>% arrange(MedRT) %>% head

# Filter criteria
with(filter_shift_sum, mean(MedRT))
with(filter_shift_sum, sd(MedRT))
with(filter_shift_sum, mean(MedRT))  + c(1.96, -1.96) * with(filter_shift_sum, sd(MedRT))

nrow(filter_shift_sum %>% filter(MedRT < 467))

goodwid_shift = as.character(filter_shift_sum %>% filter(MedRT >= 500) %>% with(WID))

shift_raw = shift_raw_inp %>% filter(WID %in% goodwid_shift)
shift_raw$LearnTest = relevel(shift_raw$LearnTest, "Weight")
levels(shift_raw$LearnTest) = c('Weight Accurate', 'Weight Inaccurate')


shift_raw_inp %>% group_by(WID, LearnTest) %>% count %>% filter(n == 180) %>% with(table(LearnTest))
shift_raw %>% group_by(WID, LearnTest) %>% count %>% filter(n == 180) %>% with(table(LearnTest))


shift_params_raw = read_json("../Modeling/anonymized_output/learn_bene/base_strats_joint_percept_params.json")
shift_strategies = data.frame(
  WID = goodwid_shift,
  Strat_SP = sapply(goodwid_shift, function(w){shift_params_raw$ind_strategies[[w]]$strat_params$sp}),
  Strat_SMP = sapply(goodwid_shift, function(w){shift_params_raw$ind_strategies[[w]]$strat_params$smp}),
  Strat_Guess = sapply(goodwid_shift, function(w){shift_params_raw$ind_strategies[[w]]$strat_params$guess})
)
rownames(shift_strategies) = NULL

shift_mod = read.csv("../Modeling/anonymized_output/learn_bene/base_strats_joint_percept_learnbene.csv")

shift_dat_all = shift_raw %>%
  #filter(IsLearning != "Learning") %>% 
  select(WID, IsLearning, Trial=TrialBase, LearnTest, Response=NormResp, TrialType, CorrectResp, Time, Prioritizes) %>% 
  merge(shift_strategies) %>% 
  merge(shift_mod %>% select(WID, Trial, ModLeft, ModBal, ModRight, LLH)) %>% 
  mutate(Response = factor(as.character(Response), levels=c('L', 'B', "R")))

shift_dat = shift_dat_all %>% filter(IsLearning != "Learning") %>% select(-IsLearning)

shift_trialdat = shift_dat %>% 
  group_by(Trial, LearnTest, TrialType) %>%
  summarize(EmpLeft = mean(Response=='L'), EmpBal = mean(Response=='B'), EmpRight = mean(Response=='R'),
            ModLeft = mean(ModLeft), ModBal = mean(ModBal), ModRight = mean(ModRight),
            NLeft = sum(Response=='L'), NBal = sum(Response=='B'), NRight = sum(Response=='R')) %>%
  mutate(EmpAcc = ifelse(TrialType=='CB', EmpBal, EmpLeft), ModAcc = ifelse(TrialType=='CB', ModBal, ModLeft)) %>%
  ungroup()

shift_strategies = shift_strategies %>% 
  merge(unique(shift_raw %>% select(WID, LearnTest))) %>% 
  mutate(WeightPct = Strat_SMP / (Strat_SP + Strat_SMP))

# In addition, there are many cases of re-used graphing and analyses which are preloaded
# elsewhere (e.g., bootstrapping correlations or graphing accuracies)
source('bb_analysis_funcs.R')

# There are also a number of bespoke functions for running the rationality analysis
source('bb_rationality_funcs.R')

# Combine all of the individual stuff into a single "core" dataset for testing vs Siegler theories
siegler_type_full = shapes_ind_full %>% subset(Shape == 'BB') %>% select(-Shape, -BaseName, -Falls) %>% mutate(Experiment = 'Shapes') %>%
  rbind(materials_ind_full %>% subset(Material == 'Pure') %>% select(-Material, -BaseName, -Falls) %>% mutate(Experiment = 'Materials')) %>%
  rbind(balance_ind_full %>% subset(Centering == 'Centered' & StrutWidth == .25 & 
                                      ((Type %in% c('Bal','CB') & Falls == 'B') | (!(Type %in% c('Bal','CB')) & Falls == 'L'))) %>% 
          select(-BaseName, -Class, -Centering, -StrutWidth, -SubClass, -Falls) %>% mutate(Experiment = 'Balance'))

###########
# NOTE: This doesn't work with CV data... need to update!!!!
#########
modpart_list = c('inc_dist', 'comb', 'add_dist', 'swap_dist', 'no_sym', 'no_weight', 'no_phys', 'just_phys', 'rules')
cv_modparts = crossval_dat %>%
  merge(crossval_dat %>% 
          filter(FitType == 'all', Strategy == 'base_strats') %>% 
          select(CVOrder, BaseFitLLH = FitLLH, BaseCVLLH = CVLLH)) %>% 
  filter(FitType == 'all', Strategy %in% modpart_list) %>% 
  mutate(dFit = FitLLH - BaseFitLLH, dCV = CVLLH - BaseCVLLH,
         Strategy = factor(Strategy, levels = modpart_list))


levels(cv_modparts$Strategy) = c('All Strategies', 'S/W/P Strategies', 'Plus Dist', 'Instead Dist', 'No Symmetry',
                              'No Weight', 'No Simulation', 'Simulation Only', 'Traditional Rules')

cv_modparts_agg = cv_modparts %>% 
  group_by(Strategy) %>% 
  summarize(AvgCV = mean(dCV), LowCV = quantile(dCV, .05), HighCV = quantile(dCV, .95),
            AvgFit = mean(dFit), LowFit = quantile(dFit, .05), HighFit = quantile(dFit, .95))

cv_vsadd1_agg = crossval_dat %>%
  mutate(IsAdd1 = substr(as.character(Strategy),1,4) == 'add1') %>% 
  merge(crossval_dat %>% 
          filter(FitType == 'all', Strategy == 'base_strats') %>% 
          select(CVOrder, BaseFitLLH = FitLLH, BaseCVLLH = CVLLH)) %>% 
  filter(FitType == 'all', IsAdd1) %>% 
  mutate(dFit = FitLLH - BaseFitLLH, dCV = CVLLH - BaseCVLLH) %>% 
  select(-IsAdd1) %>% 
  group_by(Strategy) %>% 
  summarize(AvgCV = mean(dCV), LowCV = quantile(dCV, .05), HighCV = quantile(dCV, .95),
            AvgFit = mean(dFit), LowFit = quantile(dFit, .05), HighFit = quantile(dFit, .95)) %>% 
  ungroup %>% 
  mutate(StratName = substring(as.character(Strategy), 6),
         Group = rep(c('G1', 'G2'), each = nrow(.)/2))

cv_fittype = crossval_dat %>% 
  merge(crossval_dat %>% 
          filter(FitType == 'all', Strategy == 'base_strats') %>% 
          select(CVOrder, BaseFitLLH = FitLLH, BaseCVLLH = CVLLH)) %>% 
  filter(Strategy == 'base_strats', FitType != "all") %>% 
  mutate(dFit = FitLLH - BaseFitLLH, dCV = CVLLH - BaseCVLLH,
         FitType = factor(as.character(FitType)))

cv_fittype_agg = cv_fittype %>% 
  group_by(FitType) %>% 
  summarize(AvgCV = mean(dCV), LowCV = quantile(dCV, .05), HighCV = quantile(dCV, .95),
            AvgFit = mean(dFit), LowFit = quantile(dFit, .05), HighFit = quantile(dFit, .95))
levels(cv_fittype_agg$FitType) = c("All", "Strategy", "Other")


cv_indcomp_agg = crossval_ind_vs_rules %>% 
  group_by(WID, Experiment, Rule) %>% 
  summarize(AvgDLLH = mean(dLLH), LowDLLH = quantile(dLLH, .05), HighDLLH = quantile(dLLH, .95),
            AvgDPhysLLH = mean(dPhysLLH), LowDPhysLLH = quantile(dPhysLLH, .05), HighDPhysLLH = quantile(dPhysLLH, .95)) %>% 
  ungroup




#+ ToDos -----------------------


#+ Intro ----------------------------------------



#+ Experiments --------------------------------------------

#' # Experiments 1-3
#' 
#' Introduce all three basic experiments at once (shapes / materials / balance) & describe how
#' changes in features can unveil differences in judgments
#'
#' ## Experiment 1: Shapes

#' Empirical accuracy (average accuracy = `r mean(shapes_ind$EmpAcc)`):

shapes_emp_acc_plot = plot_raw_acc(shapes_agg, 'Shape','Type',T)
print(shapes_emp_acc_plot)

#' There is a difference between how people respond to the pure blocks and the shapes. When judging the pure blocks versus shapes,
#' participants were more accurate on the balance and weight trials, but less accurate on the conflicting distance trials.
#' This suggests that participants were more likely to apply simple rules (e.g., symmetry and weight rules) with the blocks.

# shapes_acc_test = glmer(data = shapes_ind_full, family = binomial,
#                         EmpAcc ~ Type*Shape + (1 | WID) + (Type | WID) + (Shape | WID) + (Type:Shape | WID) + (1 | BaseName))

shapes_acc_test = aov(data = shapes_agg, EmpAcc ~ Type*Shape)

print(summary(shapes_acc_test))

#' ## Experiment 2: Materials
#' 
#' Empirical accuracy (average accuracy = `r mean(materials_ind$EmpAcc)`):

mat_emp_acc_plot = plot_raw_acc(materials_agg, 'Material','Type',T)
print(mat_emp_acc_plot)

#' Here there is almost no difference in accuracy between beams made of a single material and beams made of multiple
#' materials. This suggests that participants are processing the stimulus in a similar way but hints that they are 
#' incorporating material properties. 
#' This is not too surprising -- we gave them approximate weights in the beginning

mat_acc_test = aov(data = materials_agg, EmpAcc ~ Type*Material)
print(summary(mat_acc_test))

#' ## Experiment 3: Balance
#' 
#' The beam rules as stated do not account for the strut on which beam rests, assuming that it is perfectly centered
#' with negligable width. While the rules can be adapted for changes in the beam strut, these should not change people's
#' judgments -- e.g., the mass heuristic should be the same for struts of any size (so long as the blocks are not over the 
#' strut) -- so if people do change their predictions as the strut is shifted or resized, then they cannot be solely
#' relying on rules.
#' 
#' We can split the balance stimuli into the two types: ones where the configuration (along with the strut) shifts
#' from centered to uncentered, and ones where the configuration stays in the same place but the strut changes size. 
#' 
#' The first set can tell us whether people account for the weight of the beam itself, and the second set can tell us
#' whether people are accounting for the size of the strut.
#' 
#' ### Strut shifting

balance_agg_shift = subset(balance_agg, Class == 'BalAdjust')

balance_shift_test = lm(EmpAcc ~ Type*Centering, data=balance_agg_shift)
anova(balance_shift_test)

#' As a reminder, these stimuli were matched in pairs that had the same block/strut configuration but one version (centered)
#' was positioned with the strut at the center of the beam, while the other version (uncentered) was shifted to the left:
#' 
#' ![Centered Shift trial](../Scenes/Images_Balance/BalAdjust_Raw_CD_0.5_1_N_01.jpeg)
#' 
#' ![Unentered Shift trial](../Scenes/Images_Balance/BalAdjust_Adj_CD_0.5_1_N_01.jpeg)
#'
#' Note: The graphs will be cleaned up for the paper but there is always a conflict between weight and distance in these
#' stimuli -- CW means that weight falls when it is centered and distance falls when it is uncentered and CD is vice versa.
#'
#' Here we want to match on configurations and test how people's judgments change as the configuration is shifted:
#' what matters is how likely people are to say 'Left' (normalized) to the same configuration depending on whether it 
#' is centered or uncentered on the beam.
#' 
#' If people are not taking the weight of the beam into account, then they should be just as likely to say left no matter 
#' where the configuration is. However, we can see below that participants were more likely to say 'left' when the configuration
#' was in the center, which is exactly what we would expect if people were taking the weight of the beam into account:

shift_comp_dat = spread(select(balance_agg_shift, BaseName, Type, Centering, StrutWidth, EmpLeft),
                        Centering, EmpLeft)
shift_emp_acc_comp_plot = ggplot(shift_comp_dat, aes(x = Centered, y = Uncentered, color=Type, shape=StrutWidth)) +
  abl + geom_point() + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  theme_bw()
print(shift_emp_acc_comp_plot)

#' A paired t-test demonstrates that there is on average a higher propensity to say 'Left' for the centered trials:
print(with(shift_comp_dat, t.test(Centered, Uncentered, paired=T)), digits=5)

#' This suggests that people are taking the weight of the beam into account -- something that comes naturally from
#' simulation but would require ad-hoc addendums to the weight and symmetry rules to explain
#'
#' Accuracy plots:
bal_shift_emp_acc_plot = plot_raw_acc(balance_agg_shift %>% mutate(Type = factor(Type, levels=c('CD','CW'))),
             'Centering', 'Type', T) +
  scale_color_manual(values = c('orange', 'brown'))

#'
#' ### Strut size
#'
#' These stimuli were matched such that all configurations were the same, but there were four different sizes of the 
#' strut that the beam could be balanced on (2.5%, 5%, 10%, and 20% of the width of the beam). Half of the beams had 
#' struts in the center, and the other half had struts off-center

balance_agg_size = subset(balance_agg, Class == 'SizeAdjust')
balance_agg_size$NNotBal = with(balance_agg_size, NLeft+NRight)
balance_agg_size$PLgnB = with(balance_agg_size, NLeft / NNotBal)
balance_agg_size$Mod_PLgnB = with(balance_agg_size, ModLeft / (ModLeft+ModRight))

balance_size_mod = lm(EmpAcc ~ Type*factor(StrutWidth)*Centering, data=balance_agg_size)

balance_agg_size_pldat = balance_agg_size %>% group_by(StrutWidth, Centering) %>%
  summarize(PBal = mean(EmpBal), PLgnB = sum(NLeft) / sum(NNotBal), 
            NTot = sum(NLeft+NBal+NRight), NNB = sum(NNotBal),
            ModBal = mean(ModBal), ModPLgnB = mean(Mod_PLgnB)) %>%
  mutate(PBal_SE = sqrt(PBal*(1-PBal)/NTot), PLgnB_SE = sqrt(PLgnB*(1-PLgnB)/NNB))

print(anova(balance_size_mod))

#' If people are using the size of the strut in their judgments, then we would expect them to be more likely to say 
#' 'Balance' as the size increases, but their judgments of left vs right (conditioned on no balance) should not be affected.
#' 
#' And indeed balance judgments increase with greater strut sizes, but nothing else:

bal_size_balance_plot = ggplot(balance_agg_size_pldat %>% ungroup %>% 
                                 mutate(StrutWidth=factor(paste(as.numeric(as.character(StrutWidth))*10,'%',sep=''),
                                                          levels = c('2.5%','5%','10%','20%'))), 
                               aes(x=StrutWidth, y=PBal, group=Centering, 
                                   fill=Centering, ymin = PBal-PBal_SE, ymax = PBal+PBal_SE)) +
  geom_bar(stat='identity',position='dodge') + geom_linerange(position=position_dodge(.9)) + 
  xlab('Pivot Width') + theme_bw() + ylab('P(Balance)') 
bal_size_balance_plot_v2 = bal_size_balance_plot + scale_fill_grey() + 
  geom_point(aes(y=ModBal), color='red',size=4,alpha=.6,position=position_dodge(.9))
print(bal_size_balance_plot)

print(summary(aov(EmpBal ~ Centering*StrutWidth, data=balance_agg_size)))

#' While where the strut is positioned affects the propensity to say that the beam will fall left, the size does not:

bal_size_lgnb_plot = ggplot(balance_agg_size_pldat %>% ungroup %>% 
                              mutate(StrutWidth=factor(paste(as.numeric(as.character(StrutWidth))*10,'%',sep=''),
                                     levels = c('2.5%','5%','10%','20%'))), 
                               aes(x=StrutWidth, y=PLgnB, group=Centering, 
                                   fill=Centering, ymin = PLgnB-PLgnB_SE, ymax = PLgnB+PLgnB_SE)) +
  geom_bar(stat='identity',position='dodge') + geom_linerange(position=position_dodge(.9)) + 
  xlab('Pivot Width') + theme_bw() + ylab('P(Left | ~Balance)') 

bal_size_lgnb_plot_v2 = bal_size_lgnb_plot + scale_fill_grey() + 
  geom_point(aes(y=ModPLgnB), color='red',size=4,alpha=.6,position=position_dodge(.9))
print(bal_size_lgnb_plot)

print(summary(aov(PLgnB ~ Centering*StrutWidth, data=balance_agg_size, weights = NNotBal)))

#' Thus people cannot simply be using weight and symmetry rules, or their predictions would not change

#' Accuracy plot
bal_size_emp_acc_plot = plot_raw_acc(balance_agg_size %>% mutate(Type = factor(Type, levels = levels(materials_agg$Type)),
                                                                 StrutWidth = factor(paste(as.numeric(as.character(StrutWidth))*10,"%",sep=''),levels=c('2.5%','5%','10%','20%'))), 
                               'StrutWidth', 'Type', T)
print(bal_size_emp_acc_plot)


#+ Gross rule classification ----------------

#' We cannot get Siegler classifications for Pivot people if we only look at the basic trials (they only get 8 ea), and
#' many people are classified as "unclassifiable"

siegler_min_acc = 26/30

st_class_dat = siegler_type_full #%>% subset(Experiment != 'Balance')
st_class_dat$R1 = with(st_class_dat, ifelse(Type %in% c('Bal','Dist'), EmpBal, ifelse(Type == 'CD', EmpRight, EmpLeft)))
st_class_dat$R2 = with(st_class_dat, ifelse(Type == 'Bal', EmpBal, ifelse(Type == 'CD', EmpRight, EmpLeft)))
st_class_dat$R4 = with(st_class_dat, ifelse(Type %in% c('Bal','CB'), EmpBal, EmpLeft))


is_R3 = function(type, el, eb, er) {
  w_i = type == 'Weight'
  d_i = type == 'Dist'
  b_i = type == 'Bal'
  cw_i = type == 'CW'
  cd_i = type == 'CD'
  cb_i = type == 'CB'
  
  # Of the basic trials, get 10/12 right
  basic_acc = mean(c(el[w_i], el[d_i], eb[b_i]))
  if (basic_acc < 10/12) {return(F)}
  
  # Of the distance problems, get 3/4 right
  dist_acc = mean(el[d_i])
  if (dist_acc < 3/4) {return(F)}
  
  # More than 4/18 deviations from weight cues in complex
  wcue_complex = mean(c(el[cw_i], el[cb_i], er[cd_i]))
  if (wcue_complex > 14/18) {return(F)}
  
  return(T)
}

st_subj_dat = st_class_dat %>% group_by(WID,Experiment) %>% summarize(R1Acc = mean(R1), R2Acc = mean(R2), R4Acc = mean(R4), N = length(Trial),
                                                                      R3Qual = is_R3(Type, EmpLeft, EmpBal, EmpRight))

st_subj_dat$SieglerRule = with(st_subj_dat,ifelse(R1Acc > R2Acc & R1Acc > R4Acc & R1Acc > siegler_min_acc, 'Rule1',
                                                  ifelse(R2Acc > R1Acc & R2Acc > R4Acc & R2Acc > siegler_min_acc, 'Rule2',
                                                         ifelse(R4Acc > R1Acc & R4Acc > R2Acc & R4Acc > siegler_min_acc, 'Rule4',
                                                                ifelse(R3Qual, 'Rule3', 'Unclassified')))))

siegler_rule_trialcounts = with(st_class_dat, table(WID, Experiment))
siegler_rule_table = with(st_subj_dat, table(Experiment, SieglerRule))
kable(siegler_rule_table)
#xtable::xtable(siegler_rule_table)

# Test for bimodality
bimodality_test_data = st_class_dat %>% filter(Type == 'Dist') %>%
  group_by(WID, Experiment) %>% summarize(Acc = mean(EmpAcc), Correct = sum(EmpAcc), N = length(Trial))
bimodality_test_data$ExpLab = with(bimodality_test_data,
                                   factor(ifelse(Experiment=='Balance','Experiment 3: Pivot',
                                                 ifelse(Experiment=='Materials', 'Experiment 2: Materials',
                                                        'Experiment 1: Shapes')),
                                          levels = c('Experiment 1: Shapes', 'Experiment 2: Materials', 'Experiment 3: Pivot')))


bimodality_test_agg = bimodality_test_data %>% group_by(ExpLab, Correct) %>%
  summarize(N = length(WID))
bimodality_test_agg = bimodality_test_agg %>%
  merge(bimodality_test_agg %>% group_by(ExpLab) %>% summarize(Total = sum(N))) %>%
  mutate(Prop = N / Total)

make_integer_breaks = function(lims) {
  return(0:lims[2])
}
bimodality_test_graph = ggplot(bimodality_test_agg, aes(x=Correct, y=Prop)) +
  geom_bar(stat='identity') + 
  facet_grid(. ~ ExpLab, scales='free_x') +
  expand_limits(x=0) + scale_x_continuous(breaks=make_integer_breaks) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  xlab('Correct Distance Trials') + ylab('Proportion of participants') +
  theme_bw() 

bimodality_extreme_shape = with(bimodality_test_agg %>% filter(ExpLab == 'Experiment 1: Shapes') %>% 
                                  filter(Correct != max(Correct), Correct != 0),
                                1-sum(Prop))
bimodality_extreme_mat = with(bimodality_test_agg %>% filter(ExpLab == 'Experiment 2: Materials') %>% 
                                filter(Correct != max(Correct), Correct != 0),
                              1-sum(Prop))
bimodality_extreme_bal = with(bimodality_test_agg %>% filter(ExpLab == 'Experiment 3: Pivot') %>% 
                                filter(Correct != max(Correct), Correct != 0),
                              1-sum(Prop))

#' Test for bimodality -- are people basically always right or always wrong for 
#' distance trials in the simple cases? 
print(bimodality_test_graph)


#' Same for weight
bimodality_test_data_w = st_class_dat %>% filter(Type == 'Weight') %>%
  group_by(WID, Experiment) %>% summarize(Acc = mean(EmpAcc), Correct = sum(EmpAcc), N = length(Trial))
bimodality_test_data_w$ExpLab = with(bimodality_test_data,
                                     factor(ifelse(Experiment=='Balance','Experiment 3: Pivot',
                                                   ifelse(Experiment=='Materials', 'Experiment 2: Materials',
                                                          'Experiment 1: Shapes')),
                                            levels = c('Experiment 1: Shapes', 'Experiment 2: Materials', 'Experiment 3: Pivot')))


bimodality_test_agg_w = bimodality_test_data_w %>% group_by(ExpLab, Correct) %>%
  summarize(N = length(WID))
bimodality_test_agg_w = bimodality_test_agg_w %>%
  merge(bimodality_test_agg_w %>% group_by(ExpLab) %>% summarize(Total = sum(N))) %>%
  mutate(Prop = N / Total)

bimodality_test_graph_w = ggplot(bimodality_test_agg_w, aes(x=Correct, y=Prop)) +
  geom_bar(stat='identity') + 
  facet_grid(. ~ ExpLab, scales='free_x') +
  expand_limits(x=0) + scale_x_continuous(breaks=make_integer_breaks) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  xlab('Correct Weight Trials') + ylab('Proportion of participants') +
  theme_bw()

bimodality_extreme_shape_w = with(bimodality_test_agg_w %>% filter(ExpLab == 'Experiment 1: Shapes') %>% 
                                    filter(Correct != max(Correct), Correct != 0),
                                  1-sum(Prop))
bimodality_extreme_mat_w = with(bimodality_test_agg_w %>% filter(ExpLab == 'Experiment 2: Materials') %>% 
                                  filter(Correct != max(Correct), Correct != 0),
                                1-sum(Prop))
bimodality_extreme_bal_w = with(bimodality_test_agg_w %>% filter(ExpLab == 'Experiment 3: Pivot') %>% 
                                  filter(Correct != max(Correct), Correct != 0),
                                1-sum(Prop))

#' Test for bimodality -- are people basically always right or always wrong for 
#' distance trials in the simple cases? 
print(bimodality_test_graph_w)


#+ Model Fits ----------------------------------------------

#' # Model of balance judgments
#' 
#' ## Model predictions
#' 
#' ### Accuracy
#' 
#' Overall accuracy correlations:
avg_mod_acc = mean(c(shapes_agg$ModAcc, materials_agg$ModAcc, balance_agg$ModAcc))
print(avg_mod_acc)
avg_emp_acc = mean(c(shapes_agg$EmpAcc, materials_agg$EmpAcc, balance_agg$EmpAcc))
print(avg_emp_acc)

print(cor(c(shapes_agg$ModAcc, materials_agg$ModAcc, balance_agg$ModAcc),
          c(shapes_agg$EmpAcc, materials_agg$EmpAcc, balance_agg$EmpAcc)))

agg_acc_table = shapes_agg %>% select(Trial, Type, ModAcc, EmpAcc) %>% mutate(Model = 'Shapes') %>%
  rbind(materials_agg %>% select(Trial, Type, ModAcc, EmpAcc) %>% mutate(Model = 'Materials')) %>%
  rbind(balance_agg %>% select(Trial, Type, ModAcc, EmpAcc) %>% mutate(Model = 'Balance')) %>%
  group_by(Type, Model) %>% summarize(Mod = mean(ModAcc), Emp = mean(EmpAcc)) %>% ungroup %>%
  gather(From, Accuracy, Mod, Emp) %>% mutate(Col = interaction(From, Model)) %>% select(-From, -Model) %>%
  spread(Col, Accuracy)

#' 
#' Model accuracy was on average only slightly below human accuracy (model: `r avg_mod_acc`, participants: `r avg_emp_acc`). 
#' It also correlates well with the model across experiment types:
#' 

agg_acc_cor_tab = data.frame(Model = c('Shapes','Materials','Pivot-Size', 'Pivot-Shift'),
                             Correlation = c(with(shapes_agg, cor(ModAcc, EmpAcc)),
                                             with(materials_agg, cor(ModAcc, EmpAcc)),
                                             with(subset(balance_agg,Class=='SizeAdjust'), cor(ModAcc, EmpAcc)),
                                             with(subset(balance_agg,Class=='BalAdjust'), cor(ModAcc, EmpAcc))),
                             EmpAcc = c(with(shapes_agg, mean(EmpAcc)),
                                        with(materials_agg, mean(EmpAcc)),
                                        with(subset(balance_agg,Class=='SizeAdjust'), mean(EmpAcc)),
                                        with(subset(balance_agg,Class=='BalAdjust'), mean(EmpAcc))),
                             ModAcc = c(with(shapes_agg, mean(ModAcc)),
                                        with(materials_agg, mean(ModAcc)),
                                        with(subset(balance_agg,Class=='SizeAdjust'), mean(ModAcc)),
                                        with(subset(balance_agg,Class=='BalAdjust'), mean(ModAcc))))
kable(agg_acc_cor_tab)

bcvm_shapes = boot_compare_vs_model(shapes_ind_full, shapes_agg, 500)
bcvm_mat = boot_compare_vs_model(materials_ind_full, materials_agg, 500)
bcvm_pivot_size = boot_compare_vs_model(balance_ind_full %>% filter(Class=='SizeAdjust'),
                                        balance_agg %>% filter(Class=='SizeAdjust'), 500)
bcvm_pivot_shift = boot_compare_vs_model(balance_ind_full %>% filter(Class=='BalAdjust'),
                                        balance_agg %>% filter(Class=='BalAdjust'), 500)
agg_acc_split_vsmod = data.frame(Model = c('Shapes','Materials','Pivot-Size','Pivot-Shift'),
                                 NoiseCeiling = c(bcvm_shapes[1],bcvm_mat[1],bcvm_pivot_size[1],bcvm_pivot_shift[1]),
                                 ModelCor = c(bcvm_shapes[2],bcvm_mat[2],bcvm_pivot_size[2],bcvm_pivot_shift[2]))


bcvm_mse_shapes = boot_mse_vs_model(shapes_ind_full, shapes_agg, 500)
bcvm_mse_mat = boot_mse_vs_model(materials_ind_full, materials_agg, 500)
bcvm_mse_pivot_size = boot_mse_vs_model(balance_ind_full %>% filter(Class=='SizeAdjust'),
                                        balance_agg %>% filter(Class=='SizeAdjust'), 500)
bcvm_mse_pivot_shift = boot_mse_vs_model(balance_ind_full %>% filter(Class=='BalAdjust'),
                                         balance_agg %>% filter(Class=='BalAdjust'), 500)
agg_acc_mse_split_vsmod = data.frame(Model = c('Shapes','Materials','Pivot-Size','Pivot-Shift'),
                                 NoiseCeiling = c(bcvm_mse_shapes[1],bcvm_mse_mat[1],bcvm_mse_pivot_size[1],bcvm_mse_pivot_shift[1]),
                                 ModelCor = c(bcvm_mse_shapes[2],bcvm_mse_mat[2],bcvm_mse_pivot_size[2],bcvm_mse_pivot_shift[2]))

#' Can we get any better than this? Split-half correlations vs. average correlations based on half of participants:
kable(agg_acc_split_vsmod)

# print(xtable::xtable(agg_acc_split_vsmod, signif=2), print.rownames=F)

#' 
#' Note: figures to be cleaned and turned into one wide figure:
#' 

shapes_agg_acccor_plot = plot_acc_cor(shapes_agg, 'Type', 'Shape')
mat_agg_acccor_plot = plot_acc_cor(materials_agg, 'Type', 'Material')
bal_agg_accor_plot = plot_acc_cor(balance_agg, 'SubClass', 'Class')

bal_size_agg_accor_plot = plot_acc_cor(balance_agg_size %>% mutate(Type = factor(Type, levels = levels(materials_agg$Type))),
                                       'Type', 'StrutWidth') 

bal_shift_agg_accor_plot = plot_acc_cor(balance_agg_shift %>% mutate(Type = factor(Type, levels = c('CD','CW'))), 'Type','Centering') +
  scale_color_manual(values = c('orange', 'brown'))


#' Shapes
print(shapes_agg_acccor_plot)

#' Materials
print(mat_agg_acccor_plot)

#' Balance
print(bal_agg_accor_plot)

print(bal_size_agg_accor_plot)
print(bal_shift_agg_accor_plot)

#' 
#' ### Model choices
#' 
#' Beyond just fitting whether people get the beam right, we want to know whether the model explains what they choose 
#' when incorrect
#' 
#' (TODO: Are there good summary statistics to use here?)
#' 
#' Percentage of participants choosing each option (to be cleaned up):

choice_tab = data.frame(
  rbind(shapes_agg %>% select(EmpLeft, ModLeft, EmpBal, ModBal, EmpRight, ModRight) %>% sapply(mean),
        materials_agg %>% select(EmpLeft, ModLeft, EmpBal, ModBal, EmpRight, ModRight) %>% sapply(mean),
        balance_agg %>% select(EmpLeft, ModLeft, EmpBal, ModBal, EmpRight, ModRight) %>% sapply(mean)),
  row.names = c('Shapes','Materials','Balance'))

kable(choice_tab)

#' All figures to be cleaned up and placed side-by-side

#' Shapes
shapes_agg_choice_plot = plot_choices(shapes_agg, 'Shape','Type',4)
print(shapes_agg_choice_plot)

#' Materials
mat_agg_choice_plot = plot_choices(materials_agg, 'Material','Type',4)
print(mat_agg_choice_plot)

#' Balance: strut shifting
bal_agg_choices_shift_plot = plot_choices(balance_agg_shift,'Centering', 'Type', 4)
print(bal_agg_choices_shift_plot)

bal_agg_choices_shift_cdcw = plot_choices(balance_agg_shift %>% subset(Type=='CD') %>%
                                            mutate(A='a'),'Centering','A',4)
bal_agg_choices_shift_cwcd = plot_choices(balance_agg_shift %>% subset(Type=='CW') %>%
                                            mutate(A='a'),'Centering','A',4)

#' Balance: strut size adjustment
bal_agg_choices_size_plot = plot_choices(balance_agg_size,'SubClass', 'Type', 4)
print(bal_agg_choices_size_plot)



#+ Is this rules? ------------

#' # Is this rules?
#' 
#' No -- first, we can see that a model that uses a mixture of rules fits less well than the rules + simulation model

nobs_shapes = nrow(shapes_ind_full)
bic_shapes = -2*sum(shapes_agg$LLH) + log(nobs_shapes) * 14
bic_shapes_rules = -2*sum(rules_shapes$LLH) + log(nobs_shapes) * 15

nobs_mat = nrow(materials_ind_full)
bic_mat = -2*sum(materials_agg$LLH) + log(nobs_mat) *14
bic_mat_rules = -2*sum(rules_materials$LLH) + log(nobs_mat) * 15

nobs_bal = nrow(balance_ind_full)
bic_bal = -2*sum(balance_agg$LLH) + log(nobs_bal) * 14
bic_bal_rules = -2*sum(rules_balance$LLH) + log(nobs_bal) * 15

bic_mod_tot = bic_shapes + bic_mat + bic_bal
bic_rules_tot = bic_shapes_rules + bic_mat_rules + bic_bal_rules

#' BIC of model: `r bic_mod_tot`
#' 
#' BIC of rules: `r bic_rules_tot`
#' 



#'
#' Overall consistency

basic_shapes = shapes_ind_full %>% subset(Shape == 'BB') %>% merge(rules_ind_shapes %>% subset(Shape == 'BB') %>%
                                                                     mutate(RM_Left = ModLeft, RM_Bal = ModBal, RM_Right = ModRight) %>% 
                                                                     select(WID, Trial, RM_Left, RM_Bal, RM_Right)) %>% 
  merge(cv_indcomp_agg %>% select(WID, Rule))
basic_mat = materials_ind_full %>% subset(Material == 'Pure') %>% merge(rules_ind_mat %>% subset(Material == 'Pure') %>% 
                                                                          mutate(RM_Left = ModLeft, RM_Bal = ModBal, RM_Right = ModRight) %>%
                                                                          select(WID, Trial, RM_Left, RM_Bal, RM_Right)) %>% 
  merge(cv_indcomp_agg %>% select(WID, Rule))
basic_bal = balance_ind_full %>% subset(Centering == 'Centered' & StrutWidth == "0.25" & 
                                          ((Falls == 'L' & !(Type %in% c('Bal','CB'))) |
                                             (Falls == 'B' & Type %in% c('Bal','CB')))) %>%
  merge(rules_ind_bal %>% subset(Centering == 'Centered' & StrutWidth == "0.25" & 
                               ((Falls == 'L' & !(Type %in% c('Bal','CB'))) |
                                  (Falls == 'B' & Type %in% c('Bal','CB'))))  %>%
          mutate(RM_Left = ModLeft, RM_Bal = ModBal, RM_Right = ModRight) %>%
          select(WID, Trial, RM_Left, RM_Bal, RM_Right)) %>% 
  merge(cv_indcomp_agg %>% select(WID, Rule))


plot_vs_rules(basic_shapes, incl_rulemod = T, strict_off = T)
plot_vs_rules(basic_mat, incl_rulemod = T, strict_off = T)
plot_vs_rules(basic_bal, incl_rulemod = T, strict_off = T)

# For looking at the actual model predictions

ggplot(complete_params, aes(x= rules_Rule, y = strat_WProp, color=Experiment)) + 
  geom_boxplot()



#' Still checking
plot_full_v_rulesonly(shapes_ind_full, rules_ind_shapes, "Type", "Shape")



plot_full_v_rulesonly(shapes_ind_full, rules_ind_shapes, "Type", "Rule")
plot_full_v_rulesonly(materials_ind_full, rules_ind_mat, "Type", "Rule")
plot_full_v_rulesonly(balance_ind_full, rules_ind_bal, "Type", "Rule")

#' Number of individual rules
complete_params %>% with(table(Experiment,rules_Rule)) %>% kable

comp_shape_rule_llh = rules_ind_shapes %>% group_by(WID, Rule) %>% summarize(AvgRuleLLH = mean(LLH)) %>%
  merge(shapes_ind_full %>% group_by(WID) %>% summarize(AvgModLLH = mean(LLH)))
comp_mat_rule_llh = rules_ind_mat %>% group_by(WID, Rule) %>% summarize(AvgRuleLLH = mean(LLH)) %>%
  merge(materials_ind_full %>% group_by(WID) %>% summarize(AvgModLLH = mean(LLH)))
comp_bal_rule_llh = rules_ind_bal %>% group_by(WID, Rule) %>% summarize(AvgRuleLLH = mean(LLH)) %>%
  merge(balance_ind_full %>% group_by(WID) %>% summarize(AvgModLLH = mean(LLH)))

comp_rule_llh_agg = comp_shape_rule_llh %>% mutate(Experiment='Shapes') %>% rbind(comp_mat_rule_llh %>% mutate(Experiment='Materials')) %>%
  rbind(comp_bal_rule_llh %>% mutate(Experiment = 'Pivot')) %>% mutate(ModBetter = AvgModLLH > AvgRuleLLH)

#' Is the rules + sim model better? (By experiment)
comp_rule_llh_agg %>% with(table(Experiment, ModBetter)) %>% kable

#' Is the rules + sim model better (By rule)
comp_rule_llh_agg %>% with(table(Rule,ModBetter)) %>% kable

#' Test vs Siegler rules

comp_rules = st_subj_dat %>% select(WID, SieglerRule, Experiment) %>% 
  merge(complete_params %>% mutate(OurRule = rules_Rule) %>% select(WID, OurRule)) %>% 
  mutate(IsSame = SieglerRule == OurRule)
#comp_rules$IsSame = with(comp_rules, SieglerRule == OurRule)

kable(comp_rules %>% with(table(SieglerRule, OurRule)))



#' Testing cross-validated rules
ind_mod_vs_rules_comp = crossval_ind_vs_rules %>% 
  group_by(WID, Experiment, Rule) %>% 
    summarize(avgLLH = mean(dLLH), lowLLH = quantile(dLLH, .05), highLLH = quantile(dLLH, .95)) %>% 
    ungroup() %>% 
    mutate(Experiment = factor(Experiment, levels=c('sh', 'mat', 'bal')))
levels(ind_mod_vs_rules_comp$Experiment) = c('E1: Shapes', 'E2: Materials', 'E3: Pivot')
ind_mod_vs_rules_comp  = ind_mod_vs_rules_comp %>% 
  arrange(Rule, avgLLH) %>% 
  mutate(idx = 1:nrow(.))

ggplot(ind_mod_vs_rules_comp, aes(x=idx, y=avgLLH, color=Experiment, ymin=lowLLH, ymax=highLLH)) +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_pointrange() +
  theme_bw()

ind_mod_vs_phys_comp = crossval_ind_vs_rules %>% 
  group_by(WID, Experiment) %>% 
  summarize(avgLLH = mean(dPhysLLH), lowLLH = quantile(dPhysLLH, .05), highLLH = quantile(dPhysLLH, .95)) %>% 
  ungroup() %>% 
  mutate(Experiment = factor(Experiment, levels=c('sh', 'mat', 'bal')))
levels(ind_mod_vs_phys_comp$Experiment) = c('E1: Shapes', 'E2: Materials', 'E3: Pivot')
ind_mod_vs_phys_comp = ind_mod_vs_phys_comp %>% 
  arrange(avgLLH) %>% mutate(idx = 1:nrow(.))

ggplot(ind_mod_vs_phys_comp, aes(x=idx, y=avgLLH, color=Experiment, ymin=lowLLH, ymax=highLLH)) +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_pointrange() +
  theme_bw()

ruleinfo = data.frame(names = c('R1', 'R2','R3','R3a','R4'),
                      y = 28) %>%
  mutate(pltnames = ifelse(names == 'R3a','R+',as.character(names)))
ruleadds = 4*(1:nrow(ruleinfo))
names(ruleadds) = ruleinfo$names
rules_ind_diff = ind_mod_vs_rules_comp %>% arrange(Rule, avgLLH) %>% 
  select(-idx) %>% 
  mutate(plotidx = 1:nrow(.)) %>% 
  mutate(plotidx = plotidx + ruleadds[Rule]) %>% 
  merge(st_subj_dat %>% select(WID, SieglerRule))

rulebreaks = sapply(ruleinfo$names[-1], function(r) {with(subset(rules_ind_diff, Rule==r), min(plotidx)-2)})
ruleinfo$mids = sapply(ruleinfo$names, function(r) {with(subset(rules_ind_diff, Rule==r), mean(plotidx))})
ind_rule_diff_plt_base = ggplot(rules_ind_diff, aes(x = plotidx, y = avgLLH, ymin = lowLLH, ymax = highLLH, color = Experiment)) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_vline(xintercept = rulebreaks, linetype = 'solid') +
  geom_linerange() + geom_point() +
  theme_bw() + xlim(5, nrow(rules_ind_diff)+max(ruleadds)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab(expression(paste(Delta,LLH (ISR - Rules))))

ind_rule_diff_plt_lab = ind_rule_diff_plt_base +
  geom_text(aes(x=mids,y=y,label=pltnames,ymin=y,ymax=y,color=NULL), data=ruleinfo) +
  theme(legend.position = "none")

rules_ind_diff %>% with(table(Rule, avgLLH > 0)) %>% kable

table(rules_ind_diff$avgLLH > 0)
table(rules_ind_diff$lowLLH > 0)
table(rules_ind_diff$highLLH < 0)
print(binom.test(73, 98, .5), digits=4)

rules_ind_diff %>% with(table(SieglerRule, highLLH < 0))
rules_ind_diff %>% filter(highLLH < 0) %>% with(table(Rule, SieglerRule))

ind_phys_diff_plt = ggplot(ind_mod_vs_phys_comp,
                           aes(x=idx, y = avgLLH, ymin = lowLLH, ymax = highLLH, color = Experiment)) +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_linerange() + geom_point() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab(expression(paste(Delta,LLH (ISR - Simulation))))


table(ind_mod_vs_phys_comp$avgLLH > 0)
table(ind_mod_vs_phys_comp$lowLLH > 0)
table(ind_mod_vs_phys_comp$highLLH < 0)
print(binom.test(85, 98, .5), digits=4)


plot_full_v_rulesonly(shapes_ind_full %>% 
                        merge(rules_ind_diff %>% mutate(RSBetter = avgLLH >0) %>%
                                select(WID, RSBetter)), 
                      rules_ind_shapes, "Type", "RSBetter")

# Check the use of 'balance' for conflict trials (ignore QP rule)
conflict_qp_test = rbind(
  rbind(
    rules_ind_shapes %>% filter(Type %in% c('CD','CW')) %>% group_by(WID) %>% 
      summarize(Acc = mean(EmpAcc), UseBal = mean(EmpBal)) %>% mutate(Exp = "Shapes"),
    rules_ind_mat %>% filter(Type %in% c('CD','CW')) %>% group_by(WID) %>% 
      summarize(Acc = mean(EmpAcc), UseBal = mean(EmpBal)) %>% mutate(Exp = "Materials")),
  rules_ind_bal %>% filter(Type %in% c('CD','CW')) %>% group_by(WID) %>% 
    summarize(Acc = mean(EmpAcc), UseBal = mean(EmpBal)) %>% mutate(Exp = "Pivot")
)


#' We also need to show that we have reasons for not using the QP rule -- where all conflict problems should say balance. Even looking at
#' the simplest problems in each experiment, few people say balance for conflict problems:
#' 
#' Shape -- maximum use of balance across individuals:
print(max((rules_ind_shapes %>% filter(Type %in% c('CB','CD','CW'), Shape == 'BB') %>% group_by(WID) %>% summarize(MeanBal = mean(EmpBal)))$MeanBal)) # 38%

#' Materials:
print(max((rules_ind_mat %>% filter(Type %in% c('CB','CD','CW'), Material == 'Pure') %>% group_by(WID) %>% summarize(MeanBal = mean(EmpBal)))$MeanBal)) # 52%

#' Pivot:
print(max((rules_ind_bal %>% filter(Type %in% c('CB','CD','CW'), StrutWidth== .25, Centering == 'Centered') %>% group_by(WID) %>% summarize(MeanBal = mean(EmpBal)))$MeanBal)) # 47%



#+ Rationality analysis -------

#' ## Rationality analysis
#' 
#' Most of this is using function from `bb_rationality_funcs.R`
#' 


bratstrat_fl = "../Modeling/anonymized_output/rationality/basic_strats.RData"
if(file.exists(bratstrat_fl)) {
  load(bratstrat_fl)
} else {
  strats_for_all = multi_plot_phys_strats()
  strats_for_simple = multi_plot_phys_strats(data = rat_dat %>% filter(Type %in% c("Bal","Dist","Weight")))
  save(strats_for_all, strats_for_simple, file=bratstrat_fl)
}

strats_for_all_calc_agg =
  strats_for_all$data %>% 
  select(-strat_msp_smp, -strat_sp, -strat_p_only, -strat_ms_sm_s_m, -strat_guess) %>% 
  gather('Strategy', 'Prop', -sym_cost, -mass_cost, -phys_cost) %>% 
  mutate(Group = ifelse(Strategy == 'sp', 'SP',
                        ifelse(Strategy %in% c('smp', 'msp'), 'SWP/WSP',
                                  ifelse(Strategy == 'p', 'Phys Only',
                                         ifelse(Strategy == 'guessing', 'Guess', 
                                                ifelse(Strategy %in% c('sm', 'ms'), 'No Phys', 'Other')))))) %>% 
  group_by(Group, sym_cost, mass_cost, phys_cost) %>% 
  summarize(PropWin = sum(Prop) / 100)


stratgroup_plots = list()
for (stratgroup in unique(strats_for_all_calc_agg$Group)) {
  snm = gsub('/', '_', gsub(' ', '_', stratgroup))
  stratgroup_plots[[snm]] = 
    indstrat_plot(strats_for_all_calc_agg %>% filter(Group == stratgroup), stratgroup) + 
    theme(plot.title = element_text(size=18, face='bold'))
}

strats_for_present_agg =
  strats_for_all_calc_agg %>% 
  mutate(Group = factor(ifelse(Group %in% c('SP', 'SWP/WSP'), 'Used IS',
                               ifelse(Group == 'Guess', 'Guess', 'Other')),
                        levels = c('Used IS', 'Other', 'Guess'))) %>% 
  group_by(Group, sym_cost, mass_cost, phys_cost) %>% 
  summarize(PropWin = sum(PropWin)) %>% 
  ungroup

stratgroup_fullagg_plots = list()
for (stratgroup in unique(strats_for_present_agg$Group)) {
  snm = gsub('/', '_', gsub(' ', '_', stratgroup))
  stratgroup_fullagg_plots[[snm]] = 
    indstrat_plot(strats_for_present_agg %>% filter(Group == stratgroup), stratgroup) + 
    theme(plot.title = element_text(size=18, face='bold'))
}


presentation_translate = c(
  's' = 'S',
  'sp' = 'SP',
  'sm' = 'SW/WS',
  'smp' = 'SWP/WSP',
  'spm' = 'SPW',
  'm' = 'W',
  'ms' = 'SW/WS',
  'msp' = 'SWP/WSP',
  'mp' = 'WP',
  'mps' = 'WPS',
  'p' = 'P',
  'ps' = 'PS',
  'psm' = 'PSW/PWS',
  'pm' = 'PW',
  'pms' = 'PSW/PWS',
  'guessing' = 'Guess'
)
stratplot_presentation_agg = 
  strats_for_all$data %>% 
  select(-strat_msp_smp, -strat_sp, -strat_p_only, -strat_ms_sm_s_m, -strat_guess) %>% 
  gather('Strategy', 'Prop', -sym_cost, -mass_cost, -phys_cost) %>% 
  mutate(Group = factor(presentation_translate[as.character(Strategy)],
                        levels=c('SP', 'SWP/WSP', 'S', 'SW/WS', 'SPW', 'W', 'WP', 'WPS', 'P','PS',
                                 'PW', 'PSW/PWS', 'Guess'))) %>% 
  group_by(Group, sym_cost, mass_cost, phys_cost) %>% 
  summarize(PropWin = sum(Prop) / 100)


plt_present_intermediate_ratplot = 
  ggplot(data=stratplot_presentation_agg %>% filter(phys_cost == 0.12, sym_cost == 0.036, mass_cost == 0.036),
         aes(x=Group, y=PropWin)) +
  geom_bar(stat='identity') +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  ylab('% best VOC') +
  xlab("Strategy") +
  theme_bw()





dratstrat_fl = "../Modeling/anonymized_output/rationality/rate_strats.RData"
if(file.exists(dratstrat_fl)) {
  load(dratstrat_fl)
} else {
  strats_rate_for_all = multi_plot_phys_strats_ror()
  strats_rate_for_simple = multi_plot_phys_strats_ror(data = rat_dat %>% filter(Type %in% c("Bal","Dist","Weight")))
  save(strats_rate_for_all, strats_rate_for_simple, file=(dratstrat_fl))
}

ratstrat_wd_fl = "../Modeling/anonymized_output/rationality/dist_strats.RData"
if(file.exists(ratstrat_wd_fl)) {
  load(ratstrat_wd_fl)
} else {
  writeLines("Doing rationality with distance... please wait")
  ratstrat_wd_calc = multi_plot_phys_strats_wdist()
  save(ratstrat_wd_calc, file=(ratstrat_wd_fl))
}
ratstrat_wd_calc_agg = 
  ratstrat_wd_calc %>% 
  gather('Strategy', 'Prop', -sym_cost, -mass_cost, -phys_cost) %>% 
  mutate(Group = ifelse(Strategy == 'sp', 'SP',
                        ifelse(Strategy %in% c('smp', 'msp'), 'SWP/WSP',
                               ifelse(Strategy %in% c('sdp', 'dsp'), 'SDP/DSP',
                                      ifelse(Strategy %in% c('smdp', 'msdp', 'mdsp'), 'SMDP/MSDP/MDSP',
                                             ifelse(Strategy %in% c('sdmp', 'dsmp', 'dmsp'), 'SDMP/DSMP/DMSP',
                                                    ifelse(Strategy == 'p', 'Phys Only',
                                                           ifelse(Strategy == 'guessing', 'Guess', 
                                                                  ifelse(Strategy %in% c('smd', 'msd','mds'), 'SMD/MSD/MDS',
                                                                         ifelse(Strategy %in% c('sdm', 'dsm', 'dms'), 'SDM/DSM/DMS',
                                                                                'Other')))))))))) %>% 
  group_by(Group, sym_cost, mass_cost, phys_cost) %>% 
  summarize(PropWin = sum(Prop) / 100)

stratgroup_dist_plots = list()
for (stratgroup in unique(ratstrat_wd_calc_agg$Group)) {
  snm = gsub('/', '_', gsub(' ', '_', stratgroup))
  stratgroup_dist_plots[[snm]] = 
    indstrat_plot(ratstrat_wd_calc_agg %>% filter(Group == stratgroup), stratgroup) +
    theme(plot.title = element_text(size=18, face='bold'))
}


stratmix_fl = "../Modeling/anonymized_output/rationality/mixes.RData"
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

plt_rat_opt_strat = strats_for_all$plot
print(plt_rat_opt_strat)
plt_rat_opt_mix = multiplot_mix()
print(plt_rat_opt_mix)

#+ Validation experiment ----------------------------------------------

#' # Experiment 4: Generalizing to novel beams
#' 
#' If our model is explaining how people are making balance beams predictions in general, then it should generalize to 
#' new situations. Although entirely new types of beams would require explaining different uncertainty in perceptions, we can
#' introduce novel situations by combining the beam modifications from prior experiments. We should be able to explain this
#' with a zero-parameter model
#' 
#' ## Experiment
#' 
#' * 192 new beams, counterbalanced by:
#'     + Equal number of beam configurations (bal, cb, dist, weight, cd, cw)
#'     + Half pure blocks, half mix of blocks and shapes (matched trials)
#'     + Half pure material, half mixed materials
#'     + Half centered strut, half uncentered strut
#'     + Half small strut (2.5%), half large strut (10%)
#' 
#' ## Model fits
#' 
#' We can measure how well the model fits as the stimuli become more complex. If we label each stimulus by how many 
#' "changes" there are from the pure traditional stimului (where a "change" can be making the blocks into shapes, 
#' mixing materials, or placing the strut off-center), we see that the model holds up well even with multiple changes 
#' (though maybe slight fall-off when everything changes?)
#' 
#' Accuracy correlation is generally high (r = `r with(comb_agg, cor(ModAcc,EmpAcc))`)
comb_modvsemp_plot = plot_acc_cor(comb_agg, 'Type', 'NConds')
print(comb_modvsemp_plot)

#' Comparison to rules
sum(comb_agg$LLH)
sum(rules_comb$LLH)

sum(rules_comb$LLH) - sum(comb_agg$LLH)

#' Accuracies are relatively unbiased across changes; correlations decrease with more changes but still retains
#' good predictive accuracy (also showing noise ceiling -- avg. split half accuracy correlations)

bcm_comb_0 = boot_compare_vs_model(comb_agg_full %>% filter(NConds == 0), comb_agg %>% filter(NConds == 0), 500)
bcm_comb_1 = boot_compare_vs_model(comb_agg_full %>% filter(NConds == 1), comb_agg %>% filter(NConds == 1), 500)
bcm_comb_2 = boot_compare_vs_model(comb_agg_full %>% filter(NConds == 2), comb_agg %>% filter(NConds == 2), 500)
bcm_comb_3 = boot_compare_vs_model(comb_agg_full %>% filter(NConds == 3), comb_agg %>% filter(NConds == 3), 500)
bcm_comb_4 = boot_compare_vs_model(comb_agg_full %>% filter(NConds == 4), comb_agg %>% filter(NConds == 4), 500)


comb_noise_ceils = comb_noise_ceils = comb_agg %>% group_by(NConds) %>% 
   summarize(AvgEmpAcc = mean(EmpAcc), AvgModAcc = mean(ModAcc), r = cor(EmpAcc, ModAcc))
comb_noise_ceils$Noise_Ceiling = c(bcm_comb_0[1], bcm_comb_1[1], bcm_comb_2[1], bcm_comb_3[1], bcm_comb_4[1])
comb_noise_ceils$Mod_Half = c(bcm_comb_0[2], bcm_comb_1[2], bcm_comb_2[2], bcm_comb_3[2], bcm_comb_4[2])
comb_noise_ceils$Condition = c("Standard", "1 Change", "2 Changes", "3 Changes", "4 Changes")

kable(comb_noise_ceils)
print(xtable::xtable(comb_noise_ceils %>% select(Condition, Noise_Ceiling, Mod_Half), signif = 2),
      include.rownames = F)

comb_choices = plot_choices(comb_agg, 'NConds','Type',4)
print(comb_choices)

#+ Ferretti replication -------------------------------
#' # Experiment 5: explaining torque differences (replicating Ferretti et al '86)

#' A finding that has often stymied rule-based explanations of balance beams has been Ferretti et al, which finds that 
#' in the conflicting trials (and even in some of the non-conflicting trials), as the differences in the torque becomes larger,
#' children are more likely to give the correct answer. This cannot be explained through a simple rule gating system (e.g. Sielger)
#' and so has often been relegated as a curiosity in the literature (KAS - though I need to double check since I believe that some of
#' the neural network models claim they can explain this).
#' 
#' This was effectively a direct replication of Ferretti *except* we used twice as many beams (144 instead of 72), tested in 
#' adults on MTurk, and used our beam without defined struts
#' 
#' Empirically, we replicate the finding that accuracy increases as the torque difference increases as well:
ferretti_emp_plot = plot_raw_acc(ferretti_agg, 'DiffType', 'Type', T)
print(ferretti_emp_plot + xlab('Difference Level'))
ferretti_mod_plot = plot_raw_acc(ferretti_agg, 'DiffType', 'Type', F)
ferretti_rules_plot = plot_raw_acc(rules_ferretti, 'DiffType', 'Type', F)

#' Just looking at the non-balanced trials (similar to Ferretti):
ferretti_raw_test = lm(EmpAcc ~ Type*DiffType, data=ferretti_agg %>% filter(DiffType != 0))
ferretti_raw_test_nofinal = lm(EmpAcc ~ Type*DiffType, data=ferretti_agg %>% filter(DiffType %in% 1:3))
print(anova(ferretti_raw_test))
print(anova(ferretti_raw_test_nofinal))

#' We then used the same zero-parameter model to explain how we expect predictions to change over these trials, and find a 
#' good correlation between model predicted accuracies and actual accuracy across trials 
#' (r = `r with(ferretti_agg, cor(ModAcc, EmpAcc))`):
#' 
ferretti_acc_cor_plot = plot_acc_cor(ferretti_agg, 'Type','DiffType')
print(ferretti_acc_cor_plot)

with(ferretti_agg, cor(ModAcc, EmpAcc))

#' Rules not so well
ferretti_rules_cor_plot = plot_acc_cor(rules_ferretti, 'Type', 'DiffType')
print(ferretti_rules_cor_plot)

with(rules_ferretti, cor(ModAcc, EmpAcc))

sum(ferretti_agg$LLH - rules_ferretti$LLH)

#' This model could also explain the choices people made across these trials decently:
#' 
ferretti_choice_plot = plot_choices(ferretti_agg, 'DiffType', 'Type', 4)
print(ferretti_choice_plot)
ferretti_choice_plot_nonbal = plot_choices(ferretti_agg %>% filter(DiffType != 0), 'DiffType', 'Type', 4)
ferretti_choice_plot_bal = plot_choices(ferretti_agg %>% filter(DiffType == 0), 'DiffType', 'Type', 4)

ferretti_rules_choice_plot_nonbal = plot_choices(rules_ferretti %>% filter(DiffType != 0), 'DiffType','Type',4)

#' Note - it does miss the hardest distance trials. This might be an indication of a distance rule in there...

#+ Pushing use of rules around --------------

#' # Experiment 6: Shifting the use of rules
#' 
#' ## Empirical predictions


shift_dat_train = shift_raw %>% filter(IsLearning == 'Learning') %>% 
  mutate(IsCorrect = as.character(NormResp) == as.character(CorrectResp))
shift_dat_train %>% group_by(LearnTest, TrialType) %>% summarize(m = mean(IsCorrect))
shift_dat_train %>% group_by(LearnTest) %>% summarize(m = mean(IsCorrect))
shift_train_mod = glmer(data=shift_dat_train, family = binomial,
                        IsCorrect ~ TrialType*LearnTest + (1|WID) + (1|TrialBase),
                        control=glmerControl(optimizer='bobyqa', 
                                             optCtrl=list(maxfun=2e6)))
shift_train_mod_base = glmer(data=shift_dat_train, family = binomial,
                        IsCorrect ~ TrialType + (1|WID) + (1|TrialBase),
                        control=glmerControl(optimizer='bobyqa', 
                                             optCtrl=list(maxfun=2e6)))
Anova(shift_train_mod)
anova(shift_train_mod_base, shift_train_mod)

shift_dat = shift_dat %>% mutate(IsCorrect = as.character(Response) == as.character(CorrectResp))
shift_dat %>% group_by(LearnTest) %>% summarize(M = mean(IsCorrect))
shift_acc_mod = glmer(data=shift_dat, family = binomial,
                      IsCorrect ~ TrialType*LearnTest + (1|WID) + (1|Trial),
                      contrasts = list(LearnTest = contr.sum(unique(shift_dat$LearnTest)),
                                       TrialType = contr.sum(unique(shift_dat$TrialType))))
summary(shift_acc_mod)
Anova(shift_acc_mod)


shift_bysubj_pri = shift_dat %>% 
  group_by(WID, LearnTest, TrialType) %>% summarize(PriWeight = mean(Prioritizes == 'Weight')) %>% 
  ungroup

shift_bysubj_pri %>% group_by(LearnTest) %>% summarize(P = mean(PriWeight))
shift_bysubj_pri %>% group_by(TrialType, LearnTest) %>% summarize(P = mean(PriWeight))

plt_shift_priority = ggplot(shift_bysubj_pri,
                            aes(x=TrialType, y=PriWeight)) +
  geom_boxplot(aes(x=TrialType, color=LearnTest)) +
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=.3),aes(color=LearnTest, group=LearnTest)) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  ylab('% choosing side with more weight') +
  xlab("") +
  labs(color = "Training Condition")
plt_shift_priority

shift_priority_change = shift_bysubj_pri %>% 
  ungroup %>% group_by(WID, LearnTest) %>% 
  summarize(Weight = mean(PriWeight))
print(t.test(Weight~ LearnTest, data=shift_priority_change), digits=5)

mod_pri_change = glmer(data=shift_dat %>% mutate(PriWght = Prioritizes=='Weight'), family=binomial,
                      PriWght ~ LearnTest*TrialType + (1|WID) + (1|Trial) + (1|WID:TrialType),
                      contrasts = list(LearnTest = contr.sum(unique(shift_dat$LearnTest)),
                                       TrialType = contr.sum(unique(shift_dat$TrialType))))
mod_pri_change_base = glmer(data=shift_dat %>% mutate(PriWght = Prioritizes=='Weight'), family=binomial,
                            PriWght ~ TrialType + (1|WID) + (1|Trial) + (1|WID:TrialType),
                            contrasts = list(LearnTest = contr.sum(unique(shift_dat$LearnTest)),
                                             TrialType = contr.sum(unique(shift_dat$TrialType))))
anova(mod_pri_change_base, mod_pri_change)
Anova(mod_pri_change, type=2)
summary(mod_pri_change)

#' ## Model comparisons
#' 
shift_modcomp_plot = plot_choices(shift_trialdat %>% rename(Type = TrialType), 'LearnTest', 'Type', 4)

modacc_plot =ggplot(shift_trialdat, aes(x=ModAcc, y=EmpAcc, color=TrialType)) +
  geom_abline(intercept = 0, slope=1, linetype='dashed') +
  geom_point(size=2) +
  theme_bw() +
  scale_x_continuous(labels = scales::percent, limits=c(0,1)) +
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  facet_grid(.~LearnTest)
print(modacc_plot)

#' Overall correlation
with(shift_trialdat, cor(ModAcc, EmpAcc))

#' Non-weight training correlation
with(shift_trialdat %>% filter(LearnTest=='Bad Cue'), cor(ModAcc, EmpAcc))

#' Weight training correlation
with(shift_trialdat %>% filter(LearnTest=='Good Cue'), cor(ModAcc, EmpAcc))

shift_modstrats_plot = ggplot(shift_strategies, aes(x=LearnTest,y=WeightPct)) +
  geom_boxplot(aes(color=LearnTest)) +
  geom_point(position=position_jitter(width=.3),aes(color=LearnTest, group=LearnTest)) +
  xlab("Training Regime") + ylab("Proportion of SWP Strategies") +
  theme_bw() +
  theme(legend.position = 'none')
print(shift_modstrats_plot)

#' Under both training regimes, there are lower usages of the weight rule than in the standard experiments,
#' perhaps because there are not the simple conditions where there are some trials in which weight is a clear
#' and unambiguous way of providing an answer

shift_strategies %>% group_by(LearnTest) %>% summarise(AvgSWP = mean(WeightPct))

#' Nonetheless, there is a much greater usage of the weight rule by participants in the weight training condition

t.test(WeightPct ~ LearnTest, data=shift_strategies)



#+ Supplemental Analysis ----------------------------------------------------------------

#' # Supplemental Analyses
#' 

#' ## Combined experiment additional slicing
#' 
#' If we look at how choices match between the model and people mapped across all of the different conditions, it 
#' can be a bit hit or miss -- but this is because there are only two trials in each of the panels so there can 
#' be a lot of noise in the data
comb_agg %>% mutate(FullComb = interaction(Materials, Centering, StrutWidth, ShapeType)) %>% plot_choices('FullComb','Type',2)
#' 
#' Instead, we can again smear across the various dimensions and see that in general the model does a good job of
#' explaining people's choices.
#' 
#' Across shape type
plot_choices(comb_agg,'ShapeType','Type',4)
#' 
#' Across materials
plot_choices(comb_agg,'Materials','Type',4)
#' 
#' Across strut adjustments
comb_agg %>% mutate(Strut = interaction(Centering, StrutWidth)) %>% plot_choices('Strut','Type',4)
#' 
#' 

#' ## Misc
#' 
#' Individual's performance on "core" trials (close to Siegler)

#kable(siegler_type_full %>% group_by(WID, Type, Experiment) %>% summarize(Acc = mean(EmpAcc)) %>% spread(Type, Acc))

ggplot(siegler_type_full %>% group_by(WID, Type, Experiment) %>% summarize(Acc = mean(EmpAcc)),
       aes(x=Acc, fill = Experiment)) + geom_histogram(binwidth = .2) + facet_grid(.~Type)




#+ APPENDIX ------------

#' Inclusion of sub-strategies

#' We assert that (a) the symmetry heuristic, mass heuristic, and physics simulation are all needed for decisions, and (b) 
#' the process is ordered perception -> weight -> symmetry -> simulation. The ordering was justified theoretically by the complexity 
#' of processing required for each step, but we can directly test the necessity and ordering within the model
#'
#' As with the individual differences, we can cross-validate to see the predictive power of (a) models that remove one 
#' or more of the sub-parts, or (b) models that reorder the different parts
#' 
#' Need to note that here the physics model is generalized to MSPRT / multiple samples so that it doesn't necessarily
#' stop at physics (however, it always fits best with a single sample when physics is last)
#'

#' By removing parts or adding distance, we find that (a) a distance heuristic does very little, and (b) removing anything from the 
#' model makes it do worse
#' 
#' 

modparts_crossval_plot = ggplot(cv_modparts_agg, aes(x=Strategy, y=AvgCV, ymax=HighCV, ymin=LowCV)) +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_pointrange() + 
  theme_bw() +
  xlab("Integrated Strategies") +
  ylab(expression(Delta * " Crossvalidated LLH")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


modparts_crossval_present_base = ggplot(
  cv_modparts_agg,
  aes(x=Strategy, y=AvgCV, ymax=HighCV, ymin=LowCV)) +
  theme_bw() +
  xlab("Integrated Strategies") +
  ylab(expression(Delta * " Crossvalidated LLH")) +
  ylim(-850, 850) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

modparts_crossval_present_0 = 
  modparts_crossval_present_base +
  geom_pointrange(aes(y=AvgCV + 9999, ymax = HighCV + 9999, ymin = LowCV + 9999))
 
modparts_crossval_present_1 = 
  modparts_crossval_present_base +
  geom_pointrange(aes(y=AvgCV + 9999, ymax = HighCV + 9999, ymin = LowCV + 9999)) +
  geom_hline(yintercept = 0, linetype = 'dashed')

modparts_crossval_present_2 = 
  modparts_crossval_present_base + 
  geom_pointrange(
    data = cv_modparts_agg %>% 
      mutate(show = Strategy %in% c('All Strategies', 'S/W/P Strategies',
                                    'Plus Dist')),
    aes(y = AvgCV + ifelse(show, 0, 9999),
        ymin = LowCV + ifelse(show, 0, 9999),
        ymax = HighCV + ifelse(show, 0, 9999))
  ) + 
  geom_hline(yintercept = 0, linetype = 'dashed')

modparts_crossval_present_3 = 
  modparts_crossval_present_base +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = 'dashed')



modparts_reduced_present_plot_base = ggplot(
  cv_modparts_agg %>% 
    filter(Strategy %in% c("All Strategies", "Simulation Only", "Traditional Rules")) %>% 
    mutate(Strategy = plyr::revalue(Strategy, c("Traditional Rules" = "Rules Only"))),
  aes(x=Strategy, y=AvgCV, ymax=HighCV, ymin=LowCV)) +
  #geom_hline(yintercept=0, linetype='dashed') +
  #geom_pointrange() + 
  theme_bw() +
  xlab("Included Strategies") +
  ylab(expression(Delta * " Crossvalidated LLH")) +
  ylim(c(-850, 850)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

modparts_reduced_present_plot_line =
  modparts_reduced_present_plot_base + geom_hline(yintercept=0, linetype='dashed')

modparts_reduced_present_plot = 
  modparts_reduced_present_plot_line + geom_pointrange()


print(modparts_crossval_plot)

modparts_crossval_plot_out = modparts_crossval_plot + xlab(NULL) + theme(axis.text.x = element_blank())


fittype_crossval_plot = ggplot(cv_fittype_agg, aes(x=FitType, y=AvgCV, ymax=HighCV, ymin=LowCV)) +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_pointrange() + 
  theme_bw() +
  xlab("Individual-varying Parameters") +
  ylab(expression(Delta * " Crossvalidated LLH")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(fittype_crossval_plot)


add1_crossval_plot = ggplot(cv_vsadd1_agg,
                            aes(x=StratName, y=AvgCV,, ymax=HighCV, ymin=LowCV)) +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_pointrange() + 
  theme_bw() +
  facet_wrap(Group ~ ., nrow=2, scales="free_x") +
  xlab("Additional Strategy") +
  ylab(expression(Delta * " Crossvalidated LLH")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(), strip.text.x = element_blank())

print(add1_crossval_plot)

#+ Geometric analysis ---------------------------------------------------

#' # Experiment 5: Physical simulation versus configural processing
#' 
#' * How do we know people are actually using physical simulation rather than a number of configural rules?
#' * If they are using physics, they should be taking mass into account
#' * Materials experiment would provide only weak evidence because predictions can be identical across geometric/full configurations
#' * Therefore we ran a separate experiment to specifically tease apart these possibilities
#' 
#' ## Experiment
#' 
#' Discuss how it's the same as materials expt but with matched configurations
#' 
#' ## Results
#' 
#' Remember: trials are configurally matched across pure & mixed, but the results are **always** different between the two
#' 
#' Will be cleaned up for the paper, but the trial types are:
#' 
#' * CBO: unstable with pure, asymmetric balance when mixed
#' * CWCB: asymmetric balance when pure, falls to the side with more weight when mixed
#' * CWCD: asymmetric falling to the side with more distance when pure, more weight when mixed
#' * WB: symmetric balance when pure, falls to the side with more weight when mixed
#' * WW: positions are the same on both sides -- more blocks on the side that falls when pure, but more weight on the side
#' with fewer blocks when mixed
#'

geomat_agg_comp = geomat_agg %>% select(BaseName, Type, Material, EmpAcc) %>% spread(Material, EmpAcc)

#' If people are using configural processing, then we would expect participants would have the exact same response to 
#' the stimulus regardless of whether the materials are the same or different. But because the stimuli were created to 
#' have different ground truth outcomes across matched trials, this implies that the correlation of accuracies across 
#' matched trials should be r = -1. Instead we find a negligable and inconsistent correlation (r = 
#' `r with(geomat_agg_comp, cor(Pure,Mixed))`):
with(geomat_agg_comp, cor(Pure,Mixed))
with(geomat_agg_comp, cor.test(Pure,Mixed))

#geomat_emp_acc_plot = plot_raw_acc(geomat_agg, "Material", "Type", T)
#print(geomat_emp_acc_plot)

geomat_emp_comp_plot = ggplot(geomat_agg_comp, aes(x = Pure, y = Mixed, color=Type)) + 
  abl + geom_point() + theme_bw() + xlim(c(0,1)) + ylim(c(0,1))
print(geomat_emp_comp_plot)

#' 
#' ## Model fits
#' 
#' The relationship across matched trials was inconsistent -- in some cases the accuracies seem to be very anti-correlated 
#' but in others people seem to get it right in both correlations.
#' 
#' However, these differences are well explained by our 
#' model: the cases where people seem to get the configuration right in the mixed trials but wrong in the pure trials include 
#' ones where the mixed trials are conflicting weight (so the mass heuristic and simulation should both point towards the 
#' correct answer), while the pure trials are conflicting distance or asymmetric balance, both of which people typically get
#' wrong.
#' 
#' Configurations that we expect to be easy (where there is symmetric balance or both sides have the same configuration but 
#' one side has more weight) are the trials that participants typically got right across both matched trials.
#' 
#' If we try to explain this behavior with our model, we can do so with a zero-parameter model, reusing the parameter fits 
#' from the materials experiment. This explains predictions very well (r = `r with(geomat_agg, cor(ModAcc, EmpAcc))`, 
#' r(pure) = `r with(subset(geomat_agg, Material=='Pure'), cor(ModAcc, EmpAcc))`, r(mixed) = 
#' `r with(subset(geomat_agg, Material=='Mixed'), cor(ModAcc, EmpAcc))`):

geomat_agg_comp_plot = plot_acc_cor(geomat_agg, 'Type', 'Material')
print(geomat_agg_comp_plot)

#' Finally, we can see that this doesn't just describe accuracy, it picks out the distribution of choices well too:

geomat_agg_choices_plot = plot_choices(geomat_agg, 'Material', 'Type', 4)
print(geomat_agg_choices_plot)

#' 
#' With the geometric model, things aren't explained as well (r = `r with(geomat_geom, cor(ModAcc, EmpAcc))`, 
#' r(pure) = `r with(subset(geomat_geom, Material=='Pure'), cor(ModAcc, EmpAcc))`, r(mixed) = 
#' `r with(subset(geomat_geom, Material=='Mixed'), cor(ModAcc, EmpAcc))`)
#' 
geomat_geom_comp_plot = plot_acc_cor(geomat_geom, "Type", "Material")
print(geomat_geom_comp_plot)

geomat_raw = geomat_raw %>% mutate(TrNum = sapply(strsplit(as.character(TrialBase),"_"), function(x){return(x[3])}),
                                   MatchBase = paste(Type,TrNum,sep="_"))
ind_matches = geomat_raw %>% select(WID, MatchBase, Materials, NormResp) %>%
  spread(Materials, NormResp) %>% mutate(Matched = Mixed == Pure)
ind_matches_acc = ind_matches %>% group_by(WID) %>% 
  summarize(Match = mean(Matched), NMatch = sum(Matched), NDiff = sum(!Matched)) %>%
  filter((NMatch + NDiff) == 70)

pure_falls = geomat_raw %>% filter(Materials=='Pure') %>% mutate(PureFall=Falls) %>% select(MatchBase,PureFall) %>% unique
agg_matches = geomat_raw %>% merge(pure_falls) %>% group_by(MatchBase, Materials, Type) %>% 
  summarize(Acc = mean(NormResp == Falls), Left = mean(NormResp=="L"), Right = mean(NormResp=="R"),
            PFAcc = mean(NormResp == PureFall))

with(agg_matches %>% select(MatchBase,Type,Materials,Acc) %>% spread(Materials, Acc), cor(Mixed,Pure))

plot_geomat_acc_emp_match = ggplot(agg_matches %>% select(MatchBase,Type,Materials,Acc) %>% spread(Materials, Acc),
                                   aes(x=Pure, y=Mixed, color=Type)) + 
  geom_point() + theme_bw() + xlab("Accuracy (Pure)") + ylab("Accuracy (Mixed)")
print(plot_geomat_acc_emp_match)



#+ Exports ------------------------------------------

plt_ppt_notsim = 
  ggplot(shapes_agg %>%
           filter(Type == 'CD', Shape == 'BB') %>%
           mutate(N = NLeft + NBal + NRight) %>% 
           select(L = EmpLeft, B = EmpBal, R = EmpRight, N) %>% 
           gather(key='Dir', value='Pct', L, B, R) %>% 
           group_by(Dir) %>% 
           summarize(Pct = sum(Pct*N)/sum(N), TotN = sum(N)) %>% 
           mutate(SE = Pct / sqrt(TotN)) %>% 
           ungroup %>% 
           mutate(Dir = factor(Dir, levels=c('L', 'B', 'R'))),
         aes(x = Dir, y = Pct, ymin = Pct - 1.96*SE, ymax = Pct + 1.96*SE)) +
  geom_bar(stat='identity') +
  geom_linerange() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  ylab('Prop. Choosing L/B/R') +
  xlab('') +
  theme_bw()

if(export) {
  
  # Save the empirical accuracy plots
  empacctheme = theme(axis.text = element_text(size=18), axis.title = element_text(size=18))
  save_figure(shapes_emp_acc_plot + empacctheme, 
              'Shapes_Emp_Acc.pdf', 4, 4)
  save_figure(mat_emp_acc_plot + empacctheme, 'Materials_Emp_Acc.pdf', 4, 4)
  save_figure(bal_size_emp_acc_plot + empacctheme,
              'Balance_Size_Emp_Acc.pdf', 4, 4)
  save_figure(bal_shift_emp_acc_plot + empacctheme, 
              'Balance_Shift_Emp_Acc.pdf', 4, 4)
  
  save_figure(shapes_agg_acccor_plot + empacctheme, 
              'Shapes_ModVsEmp_Acc_Agg.pdf', 4, 4)
  save_figure(mat_agg_acccor_plot + empacctheme, 
              'Materials_ModVsEmp_Acc_Agg.pdf', 4, 4)
  
  save_figure(bal_size_balance_plot, 
              'Balance_Size_Bal.pdf',5,4,F)
  save_figure(bal_size_lgnb_plot, 
              'Balance_Size_LGNB.pdf',5,4,F)
  
  save_figure(bal_size_agg_accor_plot + empacctheme,
              'Balance_Size_ModVsEmp_Acc_Agg.pdf', 4, 4)
  save_figure(bal_shift_agg_accor_plot + empacctheme,
              'Balance_Shift_ModVsEmp_Acc_Agg.pdf', 4, 4)
  
  # save_figure(shapes_agg_acc_plot, 'Shapes_Agg_Acc.pdf', 4, 4)  #
  # save_figure(shapes_ind_acc_plot, 'Shapes_Ind_Acc.pdf', 4, 4)  #
  # save_figure(shapes_ind_acccor_plot, 'Shapes_ModVsEmp_Acc_Ind.pdf', 4, 4)  #
  # 
  # save_figure(mat_agg_acc_plot, 'Materials_Agg_Acc.pdf', 4, 4)  #
  # save_figure(mat_ind_acc_plot, 'Materials_Ind_Acc.pdf', 4, 4)  #
  # save_figure(mat_ind_acccor_plot, 'Materials_ModVsEmp_Acc_Ind.pdf', 4, 4)  #
  # 
  
  save_figure(shapes_agg_choice_plot, 'Shapes_Choice.pdf', 4, 6)
  save_figure(mat_agg_choice_plot, 'Materials_Choice.pdf', 2.75, 6)
  save_figure(bal_agg_choices_size_plot, 'Pivot_Size_Choice.pdf', 10, 4.16)
  save_figure(bal_agg_choices_shift_plot, 'Pivot_Shift_Choice.pdf', 2.75, 2.4)
  save_figure(bal_agg_choices_shift_cwcd, 'Pivot_Shift_Choice_Left.pdf', 2.78, 1.4)
  save_figure(bal_agg_choices_shift_cdcw, 'Pivot_Shift_Choice_Right.pdf', 2.78, 1.4)
  
  save_figure(plt_ppt_notsim, 'PPT_NotSim.pdf', 3,2)
  
  save_figure(modparts_crossval_plot, 'ModParts_CV.pdf', 5, 4)
  save_figure(fittype_crossval_plot, "FitType_CV.pdf", 3, 4)
  save_figure(add1_crossval_plot, "Add1_CV.pdf", 6, 4)
  
  #save_figure(modparts_reduced_present_plot_base, "ModParts_Red_Present_B.pdf", 3, 4)
  #save_figure(modparts_reduced_present_plot_line, "ModParts_Red_Present_L.pdf", 3, 4)
  #save_figure(modparts_reduced_present_plot, "ModParts_Red_Present_F.pdf", 3, 4)
  save_figure(modparts_crossval_present_0, "ModParts_Present_0.pdf", 5, 4)
  save_figure(modparts_crossval_present_1, "ModParts_Present_1.pdf", 5, 4)
  save_figure(modparts_crossval_present_2, "ModParts_Present_2.pdf", 5, 4)
  save_figure(modparts_crossval_present_3, "ModParts_Present_3.pdf", 5, 4)
  
  save_figure(comb_modvsemp_plot, 'Comb_ModVsEmp.pdf', 4, 4)
  save_figure(comb_choices, 'Comb_Choice.pdf', 5, 4)
  
  save_figure(ferretti_emp_plot + empacctheme, 'Ferretti_Emp_Acc.pdf', 4, 4)
  save_figure(ferretti_mod_plot + empacctheme, 'Ferretti_Mod_Acc.pdf', 4, 4)
  save_figure(ferretti_rules_plot + empacctheme, 'Ferretti_Rules_Acc.pdf', 4, 4)
  save_figure(ferretti_acc_cor_plot + empacctheme, 'Ferretti_Acc_Cor.pdf', 4, 4)
  save_figure(ferretti_acc_cor_plot, 'Ferretti_Cor_Legend.pdf', 5, 4, remove_legend=F)
  save_figure(ferretti_rules_cor_plot + empacctheme, "Ferretti_Rules_Cor.pdf", 4, 4)
  save_figure(ferretti_choice_plot_nonbal, 'Ferretti_Choice_NonBal.pdf', 5.25, 4.16)
  save_figure(ferretti_choice_plot_bal, 'Ferretti_Choice_Bal.pdf', 1.4, 2.4)
  
  # save_figure(plot_geomat_acc_emp_match, "Geomat_Acc_match.pdf", 4, 4)
  # save_figure(geomat_agg_comp_plot, "Geomat_ModCor.pdf", 4, 4)
  # save_figure(geomat_geom_comp_plot, "Geomat_GeomCor.pdf", 4, 4)
  
  save_figure(ind_rule_diff_plt_lab, "IndRuleDiff_Plot.pdf", 5, 3)
  save_figure(ind_rule_diff_plt_base, "IndRuleDiff_Legend.pdf", 6, 3, remove_legend = F)
  save_figure(ind_phys_diff_plt, "IndPhysDiff_Plot.pdf", 4, 3)
  
  save_figure(geomat_emp_comp_plot + 
                xlab('Accuracy (Pure)') + ylab('Accuracy (Mixed)'),
              "Geomat_Acc_match.pdf", 4, 4)
  save_figure(geomat_agg_comp_plot, "Geomat_ModCor.pdf", 4, 4)
  save_figure(geomat_geom_comp_plot, "Geomat_GeomCor.pdf", 4, 4)
  
  ggsave("Bimodality_Graph.png", bimodality_test_graph, width=7, height=4, units='in')
  
  save_figure(plt_rat_opt_mix + theme(axis.text = element_text(size=14), axis.title = element_text(size=14),
                                      strip.text = element_text(size=14), legend.text = element_text(size=14)), 
              "Rationality_Mixture.pdf", 10, 8, remove_legend=F)
  save_figure(plt_rat_opt_strat + theme(axis.text = element_text(size=14), axis.title = element_text(size=14),
                                        strip.text = element_text(size=14)), 
              "Rationality_Strategy.pdf", 10, 8, remove_legend=F)
  
  save_figure(plt_shift_priority, "Shift_Priority.pdf", 6, 4, remove_legend=F)
  
  save_figure(plt_present_intermediate_ratplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
              "Present_Rationality_Intermediate.pdf", 6, 4)
  
  
  sgp_theme = theme(legend.position='none', plot.title = element_text(size=24, face='bold'),
                    axis.title = element_blank())
  for (sg in names(stratgroup_plots)) {
    save_figure(stratgroup_plots[[sg]] + sgp_theme, 
                paste('RatStrat_Base_', sg, '.pdf', sep=''), 6, 6, remove_legend = F)
  }
  for (sg in names(stratgroup_dist_plots)) {
    save_figure(stratgroup_dist_plots[[sg]] + sgp_theme, 
                paste('RatStrat_Dist_', sg, '.pdf', sep=''), 6, 6, remove_legend = F)
  }
  save_figure(stratgroup_dist_plots[[1]] + labs(fill="% Best"), 'RatStrat_LEGEND.pdf', 6, 6, remove_legend=F)
  
  for (sg in names(stratgroup_fullagg_plots)) {
    save_figure(stratgroup_fullagg_plots[[sg]] + sgp_theme, 
                paste('Present_RatStrat_FullAgg_', sg, '.pdf', sep=''), 6, 6, remove_legend = F)
  }
  
  #save_figure(plot_comb_rulediff, "Comb_IndRuleDiff_Plot.pdf", 6,4)
  
}

