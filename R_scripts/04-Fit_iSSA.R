
#################################################
###                                           ###
### Fit models:                               ###
### - Step length by pre- and post- calving   ###
###   according to GC level                   ###
### - Cover and crop habitats pre- and post-  ###
###   calving according to GC levels          ###
###                                           ###
### Models are fit n times after resampling   ###
### at cluster level (cluster = individual    ###
### x gc level), and stored in                ###
### 'output/model_boots/' for later           ###
### 95% CI bootstrapping                      ###
###                                           ###
#################################################

library(sf)
library(tidyverse)

# Load data
model_dat <- readRDS('alt_output/model_dat.rds') %>%
  na.omit() %>%
  # Factor pre- and post-calving periods
  mutate(period = factor(period, 
                         levels = c('pre_calv', 'post_calv'))) %>%
  # Group for bootstraps 
  # (IMPORTANT: loading glmmTMB first will interfere with cur_group_id() 
  # function, preventing bootstrap from working. If glmmTMB is loaded, confirm 
  # that each id x cort level was assigned to a different group.)
  group_by(id, cort_ng_g) %>%
  mutate(group_id = cur_group_id()) %>%
  # Filter out elk with >3 GC samples
  filter(! id %in% c('ER_E_31', 'ER_E_20'))

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

# Create folder to store bootstrap samples
if(!dir.exists('alt_output/model_boots')) dir.create('alt_output/model_boots')

# Set number of iterations
it_n <- 1000

# Set start iteration to zero or next iteration
existing_its <- list.files('alt_output/model_boots')
iteration <- ifelse(length(existing_its) == 0, 0, 
                    max(as.numeric(str_extract_all(existing_its, '[0-9]+'))))


# Repeat loop n times to save n models
repeat {
  # Set iteration
  iteration <- iteration + 1
  # Take bootstrap sample from data
  boot_IDs <- data.frame(group_id = sample(unique(model_dat$group_id),
                                           size = length(unique(model_dat$group_id)),
                                           replace = T))
  # Build bootstrap dataframe
  model_boot <- left_join(model_dat, boot_IDs)
  # Run models
  for(i in c('crop', 'cover')) {
    model_covs <- get(paste(i, 'covs', sep = '_'))
    # Set number of random parameters for mapping
    nvar_parm <- ifelse(i %in% c('crop', 'cover'), 3, 1)
    # Set up model without fitting
    model_form <- suppressWarnings(
      glmmTMB(reformulate(model_covs, response = 'case_'),
              family = poisson(), 
              map = list(theta = factor(c(NA, 1:nvar_parm))),
              data = model_boot, doFit = F))
    # Set variance of random intercept to 10^6
    model_form$parameters$theta[1] <- log(1e6)
    # Fit model using large fixed variance
    model_fit <- glmmTMB:::fitTMB(model_form)
    
    # Save models
    saveRDS(model_fit, paste('alt_output/model_boots/', 
                             paste(i, 'model', iteration, sep = '_'), 
                             '.rds', sep =''))
  }
  # When iterations reach n, break the script
  if(iteration == it_n) break
  
}
