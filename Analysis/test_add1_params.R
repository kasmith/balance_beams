library(tidyverse)
library(ggstance)

dat_all = read.csv('../Modeling/output/strat_crossval/all_crossval.csv')
dat_add1 = read.csv('../Modeling/output/strat_crossval/add1_crossval.csv')

dat_full = dat_add1 %>% merge(
  dat_all %>% filter(Strategy == 'base_strats', FitType == 'all') %>% 
    select(CVOrder, FitBase = FitLLH, CVBase = CVLLH)
) %>% mutate(
  dLLH_Base = FitBase - FitLLH, dLLH_CV = CVBase - CVLLH
)

dat_agg = dat_full %>% 
  group_by(Strategy) %>% 
  summarize(dLLHBase_avg = mean(dLLH_Base), dLLHBase_low = quantile(dLLH_Base, .1), dLLHBase_high = quantile(dLLH_Base, .9),
            dLLHCV_avg = mean(dLLH_CV), dLLHCV_low = quantile(dLLH_CV, .1), dLLHCV_high = quantile(dLLH_CV, .9),
            Strat_avg = mean(AddStratParam), Strat_low = quantile(AddStratParam, .1), Strat_high = quantile(AddStratParam, .9))


ggplot(data = dat_agg,
       aes(x = Strat_avg, y=dLLHCV_avg)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_linerangeh(aes(xmin = Strat_low, xmax= Strat_high), color='grey') +
  geom_linerange(aes(ymin = dLLHCV_low, ymax = dLLHCV_high), color='grey') +
  geom_point() +
  theme_bw()

dat_agg %>% filter(Strat_avg > .1) %>% arrange(-dLLHCV_avg)
