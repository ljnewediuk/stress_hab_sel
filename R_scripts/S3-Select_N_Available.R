#################################################
###                                           ###
### Test coefficient estimates under various  ###
### numbers of available:used point ratios to ###
### determine minimum number to estimate      ###
### coefficients                              ###
###                                           ###
#################################################

library(tidyverse)
library(sf)
library(raster)
library(amt)
library(glmmTMB)

# Load data
dat <- readRDS('output/stress_dat.rds')
cover <- raster::raster('rasters/cover_rast.tif')
names(cover) <- 'cover'

# Make track
trk <- amt::mk_track(dat, 
                     .x=X, .y=Y, .t=dat_time, id=animal_ID, 
                     # Keep cort and stress response data
                     cort_ng_g, period, stress_resp,
                     crs = sp::CRS("+init=epsg:26914"))

# Extract random steps (1x, 5x, 10x, 20x, 30x, 40x 50x, 100x used steps)
avail_samples <- data.frame()

if(file.exists('output/N_avail_df.rds')) {
  
  avail_samples <- readRDS('output/N_avail_df.rds')
  
} else {
  
  for(n_steps in c(1, 5, 10, 20, 30, 40, 50, 100, 500, 1000)) {
    
    stps <- track_resample(trk, rate=hours(2), tolerance=minutes(30)) %>%
      # Make sure bursts have at least three points
      filter_min_n_burst(min_n = 3) %>% 
      steps_by_burst(keep_cols = 'start') %>% 
      random_steps(n = n_steps) %>%
      # Extract land cover
      extract_covariates(cover, where = 'end') %>%
      group_by(id) %>%
      # Scale cort by ID
      mutate(log_sl_ = log(sl_ + 1),
             cos_ta_ = cos(ta_),
             cort_ng_g_sc = scale(cort_ng_g, scale = T, center = F),
             period = factor(period, 
                             levels = c('pre_calv', 'post_calv'))) %>%
      na.omit() %>%
      group_by(id, cort_ng_g) %>%
      mutate(group_id = cur_group_id()) %>%
      # Filter out elk with >3 GC samples
      filter(! id %in% c('ER_E_31', 'ER_E_20'))
    
    # Set max iterations to 10^15 to aid convergence
    glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))
    
    # Cover model
    cover_covs <- c('cort_ng_g_sc:cover:period', 'cort_ng_g_sc:cover',
                    'cover:period', 'cover', 'log_sl_', 'cos_ta_',
                    '(1 | step_id_)', '(0 + cort_ng_g_sc + cover | id)')
    
    # Set up model without fitting
    model_form <- suppressWarnings(
      glmmTMB(reformulate(cover_covs, response = 'case_'),
              family = poisson(), 
              map = list(theta = factor(c(NA, 1:3))),
              data = stps, doFit = F))
    # Set variance of random intercept to 10^6
    model_form$parameters$theta[1] <- log(1e6)
    # Fit model using large fixed variance
    model_fit <- glmmTMB:::fitTMB(model_form) %>%
      broom.mixed::tidy() %>%
      mutate(n_steps) %>%
      filter(term %in% c('cover', 'log_sl_', 'cos_ta_', 
                         'cover:cort_ng_g_sc','cover:periodpost_calv', 
                         'cover:cort_ng_g_sc:periodpost_calv')) %>%
      select(n_steps, term, estimate, std.error)
    
    avail_samples <- rbind(avail_samples, model_fit)
    
  }
  
  # Save
  saveRDS(avail_samples, 'output/N_avail_df.rds')
  
}

# Plot coefficient estimates by used:available points
ggplot() +
  geom_pointrange(data = avail_samples, aes(x = factor(n_steps), y = estimate,
                                            ymin = estimate - std.error*2,
                                            ymax = estimate + std.error*2)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_wrap(~ term, scales = 'free')

avail_samples <- rbind(avail_samples, avail_samples2)


