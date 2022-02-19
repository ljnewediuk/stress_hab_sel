
#################################################
###                                           ###
### Validate models using UHC plots           ###
###                                           ###
### Fit model to only individuals with both   ###
### pre- and post-calving data                ###
###                                           ###
#################################################

library(sf)
library(tidyverse)
library(survival)

# Load data
dat <- readRDS('output/model_dat.rds') %>%
  na.omit() 

cover <- c('cover', 'cover:cort_ng_g_sc', 'log_sl_', 'cos_ta_')
crop <-c('crop', 'crop:cort_ng_g_sc', 'log_sl_', 'cos_ta_')

model_summs <- data.frame()

for(i in c(15, 16, 19, 21, 23, 28)) {
  for(per in c('pre_calv', 'post_calv')) {
    for(covs in c('cover', 'crop')) {
      
      indiv <- paste0('ER_E_', i)
      
      indiv_dat <- dat %>%
        filter(id == indiv & period == per)
      
      # Fit iSSA
      try({
        indiv_mod <- clogit(reformulate(c(get(covs), 'strata(step_id_)'), 'case_'), 
                            data= indiv_dat)
      })
      
      # Tidy
      tidy_mod <- broom::tidy(indiv_mod) %>%
        mutate(indiv, covs, per) %>%
        select(indiv, covs, per, term, estimate, std.error)
      
      if(any(tidy_mod$std.error > 50)) next
      
      model_summs <- rbind(model_summs, tidy_mod)  
      
    }
  }
}

sub <- model_summs %>%
  filter(covs == 'cover' & per == 'pre_calv')

ggplot(data = model_summs, aes(x = term, y = estimate, group = indiv)) +
  geom_point() +
  # geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) +
  facet_grid(rows = vars(covs), cols = vars(per))
  











