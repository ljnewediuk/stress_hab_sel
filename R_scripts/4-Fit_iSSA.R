
#################################################
###                                           ###
### Fit models:                               ###
### - Step length by pre- and post- calving   ###
###   according to GC level                   ###
### - Cover and crop habitats pre- and post-  ###
###   calving according to GC levels          ###
###                                           ###
#################################################

library(sf)
library(tidyverse)
library(amt)
library(glmmTMB)

# Load data
model_dat <- readRDS('output/model_dat.rds') %>%
  na.omit()

# Set max iterations to 10^15 to aid convergence
glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))

# Define covariates for models
# Step length model
sl_covs <- c('log_sl_:cort_ng_g_sc:period', 'log_sl_:period', 'log_sl_',
             '(1 | step_id_)', '(0 + cort_ng_g_sc | id)')
# Cover model
cover_covs <- c('cort_ng_g_sc:cover:period', 'cort_ng_g_sc:cover', 
              'cover:period', 'cover',
              '(1 | step_id_)', '(0 + cort_ng_g_sc + cover | id)')
# Crop model
crop_covs <- c('cort_ng_g_sc:crop:period', 'cort_ng_g_sc:crop', 
                'crop:period', 'crop',
                '(1 | step_id_)', '(0 + cort_ng_g_sc + crop | id)')

hab_covs <- c('crop', 'cover', 'crop:period', 'cover:period',
              '(1 | step_id_)', '(0 + cover + crop | id)')

# Run models
for(i in c('crop', 'cover', 'sl', 'hab')) {
  model_covs <- get(paste(i, 'covs', sep = '_'))
  # Set up model without fitting
  model_form <- glmmTMB(reformulate(model_covs, response = 'case_'),
                      family=poisson(), 
                      data = model_dat, doFit=FALSE)
  # Set variance of random intercept to 10^6
  model_form$parameters$theta[1] <- log(1e6)
  nvar_parm <- length(model_form$parameters$theta)
  model_form$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
  # Fit model using large fixed variance
  model_fit <- glmmTMB:::fitTMB(model_form)
  # Assign
  assign(paste(i, 'model', sep = '_'), model_fit)
}

# Plot step length by cort
stps_summ <- model_dat %>%
  as.data.frame() %>%
  ungroup() %>%
  filter(case_ == TRUE) %>%
  # Summarize mean step length at each cort value by period
  group_by(id, cort_ng_g, period) %>%
  summarize(mean_sl = mean(sl_, na.rm = T))
# Plot
ggplot(stps_summ, aes(x = cort_ng_g, y = log(mean_sl), col = period)) + 
  geom_point() +
  geom_smooth(method = 'lm')

# Plot habitat models
for(i in c('crop', 'cover', 'hab')) {
  tidy_hab_model <- broom.mixed::tidy(get(paste(i, 'model', sep = '_'))) %>%
    # Get confidence intervals
    mutate(lower = estimate - std.error*1.96,
           upper = estimate + std.error*1.96) %>%
    # Filter out unneeded coefficients
    filter(! term == '(Intercept)') %>%
    filter(! effect == 'ran_pars')
  # Plot
  hab_plot <- ggplot(tidy_hab_model, aes(x = term, y = estimate)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper)) + 
    theme(axis.text.x = element_text(angle = 45))
  # Assign
  assign(paste(i, 'plot', sep = '_'), hab_plot)
}


# Save model outputs
for(i in c('crop', 'cover', 'sl', 'hab')) {
  saveRDS(get(paste(i, 'model', sep = '_')), 
          paste('output/', i, '_model_results.rds', sep = ''))
}


