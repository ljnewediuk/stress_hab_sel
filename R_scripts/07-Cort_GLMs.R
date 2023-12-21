
#################################################
###                                           ###
### Fit brms models to test whether GC levels ###
### change after calving (GC ~  period)       ###      
###                                           ###
#################################################

library(tidyverse)
library(brms)
library(emmeans)
library(tidybayes)
library(performance)
library(bayestestR)

# Load hormone data with calving dates
dat <- readRDS('output/horm_calv_dat.rds') %>%
  mutate(cort_ng_g_sc = scale(cort_ng_g)[,1])

# Fit model to test if glucocorticoids higher or lower after calving
mod <- brm(cort_ng_g_sc ~ period + factor(year) + (1 | animal_ID),
           data = dat, family = gaussian, 
           iter = 10000, warmup = 5000, chains = 4, cores = 4, 
           prior = prior(normal(0,1), class = b),
           control = list(adapt_delta = 0.99, max_treedepth = 18),
           backend = 'cmdstanr')
  
# Get mean and sd of cort for unscaling
mean_cort <- mean(dat$cort_ng_g, na.rm = T)
sd_cort <- sd(dat$cort_ng_g, na.rm = T)

# Get draws of marginal estimates
emmeans_mod <- emmeans(mod, specs = 'period')
draws_mod <- gather_emmeans_draws(emmeans_mod, value = 'estimate') %>%
  mutate(estimate_unsc = estimate * sd_cort + mean_cort)

# Save draws for plotting
saveRDS(draws_mod, 'output/cort_model_draws.rds')

