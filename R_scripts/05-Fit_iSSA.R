
#################################################
###                                           ###
### Fit models for distance to cover based on ###
### cort levels, period with respect to       ###
### parturition date, and days since calving  ###
###                                           ###
#################################################

library(tidyverse)

# Load data
model_dat <- readRDS('derived_data/issa_model_dat.rds') %>%
  # na.omit() %>%
  # Scale distance data
  mutate(across(c(dist_to_crop, dist_to_cover, cort_ng_g), function(x) as.vector(scale(x))),
         # Add log(SL) and cos(TA)
         log_sl_ = log(sl_ + 1),
         cos_ta_ = cos(ta_))

# Sample bootstrap data (leave-one-individual-out) for CIs around RSS
bs_list <- list()
for(i in 1:length(unique(model_dat$id))) {
  
  bs_samp <- model_dat %>%
    filter(! id == unique(model_dat$id)[i])
  bs_list[[i]] <- bs_samp
  
}

library(glmmTMB)

# Set max iterations to 10^15 to aid convergence
glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))

# Basic model covariates
model_covs <-  c('dist_to_cover',
                 'dist_to_cover:cort_ng_g',
                 '(1 | step_id_)',
                 '(0 + dist_to_cover | id)',
                 '(0 + dist_to_cover:cort_ng_g |id)')

# Add covs for "days since calving" model
model_covs_d <-  c(model_covs,
                'dist_to_cover:cort_ng_g:d_to_calv',
                '(0 + dist_to_cover:cort_ng_g:d_to_calv |id)')
# For "period" model
model_covs_p <-  c(model_covs,
                  'dist_to_cover:cort_ng_g:period')
# Full model
model_covs_f <-  c(model_covs,
                   'dist_to_cover:cort_ng_g:d_to_calv',
                   '(0 + dist_to_cover:cort_ng_g:d_to_calv |id)',
                   'dist_to_cover:cort_ng_g:period')

# Function to fit models
fit_mod <- function(covs, nvar_parm, dat) {
  # Set up model without fitting
  model_form <- suppressWarnings(
    glmmTMB(reformulate(covs, response = 'case_'),
            family = poisson(), 
            map = list(theta = factor(c(NA, 1:nvar_parm))),
            data = dat, doFit = F))
  # Set variance of random intercept to 10^6
  model_form$parameters$theta[1] <- log(1e6)
  # Fit model using large fixed variance
  model_fit <- glmmTMB:::fitTMB(model_form)
  # Return the glmmTMB object
  return(model_fit)
}

# Fit "full" model
model_f <- fit_mod(model_covs_f, 3, model_dat)

# Fit bootstrap models in list
bs_mods <- list()
for(i in 1:length(bs_list)) {
  
  dat <- bs_list[[i]]
  model_f <- fit_mod(model_covs_f, 3, dat)
  bs_mods[[i]] <- model_f
  
}


# Save model output (mean and bootstraps)
saveRDS(model_f, 'output/model_fit_f.rds')
saveRDS(bs_mods, 'output/model_fit_bs.rds')

