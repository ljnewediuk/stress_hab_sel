
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

library(tidyverse)

# Load data
model_dat <- readRDS('derived_data/issa_model_dat.rds') %>%
  # na.omit() %>%
  # Scale distance data
  mutate(across(c(dist_to_crop, dist_to_cover, cort_ng_g), function(x) as.vector(scale(x))),
         # Make the post-calving period only ~ 30 d when calf most vulnerable
         period = ifelse(d_to_calv > 3 | d_to_calv < -16, 'pre', 'post'),
         # Factor land cover
         across(c(forest, crop, cover), factor),
         # Add log(SL) and cos(TA)
         log_sl_ = log(sl_ + 1),
         cos_ta_ = cos(ta_))

library(glmmTMB)

# Set max iterations to 10^15 to aid convergence
glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))

# Full model covariates
model_covs <-  c('I(log_sl_)',
                'dist_to_cover',
                'dist_to_cover:cort_ng_g',
                'dist_to_cover:cort_ng_g:period',
                '(1 | step_id_)',
                '(0 + I(log_sl_) | id)',
                '(0 + dist_to_cover | id)',
                '(0 + dist_to_cover:cort_ng_g |id)')

# Set number of random parameters for mapping
nvar_parm <- 3
# Set up model without fitting
model_form <- suppressWarnings(
  glmmTMB(reformulate(model_covs, response = 'case_'),
          family = poisson(), 
          map = list(theta = factor(c(NA, 1:nvar_parm))),
          data = model_dat, doFit = F))
# Set variance of random intercept to 10^6
model_form$parameters$theta[1] <- log(1e6)
# Fit model using large fixed variance
model_fit <- glmmTMB:::fitTMB(model_form)

# Save model output
saveRDS(model_fit, 'output/model_fit.rds')








# Create folder to store bootstrap samples
if(!dir.exists('output/model_boots')) dir.create('output/model_boots')

# Set number of iterations
it_n <- 500

# Set start iteration to zero or next iteration
existing_its <- list.files('output/model_boots')
iteration <- ifelse(length(existing_its) == 0, 0, 
                    max(as.numeric(str_extract_all(existing_its, '[0-9]+'))))


# Repeat loop n times to save n models
repeat {
  # Set iteration
  iteration <- iteration + 1
  # Take bootstrap sample from data
  boot_IDs <- sample(unique(model_dat$uid), 
                     size = (length(unique(model_dat$uid))-1), replace = F)
  # Build bootstrap dataframe
  model_boot <- model_dat %>%
    filter(uid %in% boot_IDs)
  # Run models
  for(i in c('full')) {
    model_covs <- get(paste(i, 'covs', sep = '_'))
    # Set number of random parameters for mapping
    nvar_parm <- 5
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
    #Skip to next iteration if convergence issues
    if (inherits(model_fit, "warning")) cat('Model convergence issue; skipping')
    if (inherits(model_fit, "warning")) next
    # Save models
    saveRDS(model_fit, paste('output/model_boots/', 
                             paste(i, 'model', iteration, sep = '_'), 
                             '.rds', sep =''))
  }
  # When iterations reach n, break the script
  if(iteration == it_n) break
  
}

