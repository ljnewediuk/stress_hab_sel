
library(tidyverse)
library(brms)
library(sf)

# Load data linking uids to labels
uid_IDs <- readRDS('output/stress_dat.rds') %>%
  ungroup() %>%
  st_drop_geometry() %>%
  select(uid, label) %>%
  distinct()

# Get assigned IDs for fecal samples
dna_IDs <- readRDS('input/final_sample_IDs.rds') %>%
  # Join uids for model
  left_join(uid_IDs) %>%
  # pull out only DNA-identified
  filter(identification_type == 'dna' & ! is.na(uid))

# Get uIDs for samples at least 5 days outside estimated calving events
out_calv_IDs <- readRDS('output/horm_calv_dat.rds') %>%
  left_join(uid_IDs) %>%
  filter(d_to_calv > 5 | d_to_calv < -5)

# Load data for iSSA
model_dat <- readRDS('derived_data/issa_model_dat.rds') %>%
  # na.omit() %>%
  # Scale distance data
  mutate(across(c(dist_to_crop, dist_to_cover, cort_ng_g), function(x) as.vector(scale(x))),
         # Add log(SL) and cos(TA)
         log_sl_ = log(sl_ + 1),
         cos_ta_ = cos(ta_)) %>%
  # Add year
  mutate(year = factor(year(t1_)))

# exclude samples assigned by model for DNA-only model
model_dat_dna <- model_dat %>%
  filter(uid %in% dna_IDs$uid)

# exclude samples within 5 d of calving
model_dat_out_calv <- model_dat %>%
  filter(uid %in% out_calv_IDs$uid)

# Load hormone data with calving dates and filter out non-gen identified
horm_dat <- readRDS('output/horm_calv_dat.rds') %>%
  mutate(cort_ng_g_sc = scale(cort_ng_g)[,1]) %>%
  filter(label %in% dna_IDs$label)

# Model covariates
model_covs <-  c('dist_to_cover',
                 'dist_to_cover:cort_ng_g',
                 '(1 | step_id_)',
                 '(0 + dist_to_cover | id)',
                 '(0 + dist_to_cover:cort_ng_g |id)',
                 'dist_to_cover:cort_ng_g:d_to_calv',
                 '(0 + dist_to_cover:cort_ng_g:d_to_calv |id)',
                 'dist_to_cover:cort_ng_g:period')
# With year added
model_covs_y <- c('year',
                  model_covs)

library(glmmTMB)

# Set max iterations to 10^15 to aid convergence
glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))

# Set DNA model without fitting
model_form_dna <- suppressWarnings(
  glmmTMB(reformulate(model_covs, response = 'case_'),
          family = poisson(), 
          map = list(theta = factor(c(NA, 1:3))),
          data = model_dat_dna, doFit = F))
# Set variance of random intercept to 10^6
model_form_dna$parameters$theta[1] <- log(1e6)
# Fit model using large fixed variance
model_fit_dna <- glmmTMB:::fitTMB(model_form_dna)

# Set out-of-calving-event-only model without fitting
model_form_out_of_calv <- suppressWarnings(
  glmmTMB(reformulate(model_covs, response = 'case_'),
          family = poisson(), 
          map = list(theta = factor(c(NA, 1:3))),
          data = model_dat_out_calv, doFit = F))
# Set variance of random intercept to 10^6
model_form_out_of_calv$parameters$theta[1] <- log(1e6)
# Fit model using large fixed variance
model_fit_out_of_calv <- glmmTMB:::fitTMB(model_form_out_of_calv)

# Set year model without fitting
model_form_y <- suppressWarnings(
  glmmTMB(reformulate(model_covs_y, response = 'case_'),
          family = poisson(), 
          map = list(theta = factor(c(NA, 1:3))),
          data = model_dat, doFit = F))
# Set variance of random intercept to 10^6
model_form_y$parameters$theta[1] <- log(1e6)
# Fit model using large fixed variance
model_fit_y <- glmmTMB:::fitTMB(model_form_y)

# Fit model to test if glucocorticoids higher or lower after calving
mod <- brm(cort_ng_g_sc ~ period + factor(year) + (1 | animal_ID),
           data = horm_dat, family = gaussian, 
           iter = 10000, warmup = 5000, chains = 4, cores = 4, 
           prior = prior(normal(0,1), class = b),
           control = list(adapt_delta = 0.99, max_treedepth = 18),
           backend = 'cmdstanr')

# Fit model to test if glucocorticoids higher or lower after calving
mod <- brm(cort_ng_g_sc ~ period + factor(year) + (1 | animal_ID),
           data = horm_dat, family = gaussian, 
           iter = 10000, warmup = 5000, chains = 4, cores = 4, 
           prior = prior(normal(0,1), class = b),
           control = list(adapt_delta = 0.99, max_treedepth = 18),
           backend = 'cmdstanr')

# Summarize
# Direction of effects are the same as the full model with both DNA-identified
# and model identified samples, suggesting the model  is robust to misidentifications.
summary(model_fit)

