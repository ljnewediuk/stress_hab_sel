
#################################################
###                                           ###
### Fit Muff model to only individuals with   ###
### data in both periods to ensure those      ###
### individuals are not driving the pattern   ###
###                                           ###
#################################################

library(sf)
library(tidyverse)
library(survival)

# Load data
dat <- readRDS('output/model_dat.rds') %>%
  na.omit() %>%
  # Filter only individuals with data in both periods
  filter(id %in% c('ER_E_15', 'ER_E_16', 'ER_E_19', 
                   'ER_E_21', 'ER_E_23', 'ER_E_28')) %>%
  mutate(period = factor(period, levels = c('pre_calv', 'post_calv'),
                         ordered = T))

# Load glmmTMB
library(glmmTMB)

# Set max iterations to 10^15 to aid convergence
glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))

# Define covariates for models
# Cover model
cover_covs <- c('cort_ng_g_sc:cover:period', 'cort_ng_g_sc:cover',
                'cover:period', 'cover', 'log_sl_', 'cos_ta_',
                '(1 | step_id_)', '(0 + cort_ng_g_sc + cover | id)')
# Crop model
crop_covs <- c('cort_ng_g_sc:crop:period', 'cort_ng_g_sc:crop', 
               'crop:period', 'crop', 'log_sl_', 'cos_ta_',
               '(1 | step_id_)', '(0 + cort_ng_g_sc + crop | id)')

for(i in c('crop', 'cover')) {
  
  # Set up model without fitting
model_form <- suppressWarnings(
  glmmTMB(reformulate(get(paste0(i, '_covs')), response = 'case_'),
          family = poisson(), 
          map = list(theta = factor(c(NA, 1:3))),
          data = dat, doFit = F))
# Set variance of random intercept to 10^6
model_form$parameters$theta[1] <- log(1e6)
# Fit model using large fixed variance
model_fit <- glmmTMB:::fitTMB(model_form)

assign(paste0(i, '_model_out'), 
       broom.mixed::tidy(model_fit) %>%
         filter(term %in% c('crop', 'cover', 'crop:cort_ng_g_sc', 'cover:cort_ng_g_sc',
                            'crop:period.L',  'crop:cort_ng_g_sc:period.L',
                            'cover:period.L',  'cover:cort_ng_g_sc:period.L',
                            'log_sl_', 'cos_ta_')) %>%
         mutate(estimate = round(estimate, digits = 2),
                CI_95 = paste0('(', round(estimate - std.error*1.96, digits = 2), ', ', 
                               round(estimate + std.error*1.96, digits = 2), ')')) %>%
         select(term, estimate, CI_95)) 

}

