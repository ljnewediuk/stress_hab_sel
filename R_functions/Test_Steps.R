#################################################
###                                           ###
### Two functions to first sample a raster at ###
### a specified number of available steps     ###
### from location data, then fit a model      ### 
### using the data.                           ### 
###                                           ###
### First function 'test_steps' returns steps ###
### data for iSSA, second function 'fit_mod'  ### 
### fits model to the steps data and returns  ###
### a data.frame with the number of steps and ###
### coefficients.                             ###
###                                           ###
#################################################

library(tidyverse)
library(sf)
library(stars)
library(amt)
library(glmmTMB)

# Test steps function
test_steps <- function(dat, n_steps) {
  
  # Make track and resample for iSSA
  steps_dat <- dat %>%
    # Add projected coordinates then drop geometry
    mutate(long = st_coordinates(dat)[, 1], lat = st_coordinates(dat)[, 2]) %>%
    st_drop_geometry() %>%
    amt::mk_track(.x = long, .y = lat, .t = dat_time, id = animal_ID, 
                  # Keep data columns
                  uid, yr, cort_ng_g, t3_ng_g, period, d_to_calv,
                  crs = sp::CRS("+init=epsg:26914")) %>%
    # Extract random steps 
    track_resample(rate=minutes(30), tolerance=minutes(5)) %>%
    # Make sure bursts have at least three points
    filter_min_n_burst(min_n = 3) %>% 
    steps_by_burst(keep_cols = 'start') %>%
    random_steps(n = n_steps)
  
  issa_dat <- data.frame()
  
  for(year in c(2019, 2020)) {
    # Subset rows containing time = either 2019 or 2020
    sub_dat <- steps_dat %>%
      filter(yr == year)
    
    dist_rast <- raster::raster(paste0('rasters/coverdistance_to_rast_', year, '.tif'))
    # Rename to cover type
    names(dist_rast)[1] <- paste0('dist_to_cover')
    # Extract land cover
    lc_vals <- sub_dat %>%
      extract_covariates(dist_rast, where = 'end') 
    # Combine data
    sub_dat <- sub_dat %>% 
      left_join(lc_vals) %>%
      # Make the post-calving period only ~ 30 d when calf most vulnerable
      mutate(across(c(dist_to_cover, cort_ng_g), function(x) as.vector(scale(x))),
             period = ifelse(d_to_calv > 3 | d_to_calv < -16, 'pre', 'post'),
             # Add log(SL)
             log_sl_ = log(sl_ + 1))
    
    issa_dat <- rbind(issa_dat, sub_dat)
    
  }
  
  return(issa_dat)
  
} 

# Fit model function
fit_mod <- function(dat) {
  # Set max iterations to 10^15 to aid convergence
  glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))
  
  # Cover model
  # Full model covariates
  model_covs <-  c('dist_to_cover',
                   'dist_to_cover:cort_ng_g',
                   'dist_to_cover:cort_ng_g:period',
                   'dist_to_cover:cort_ng_g:d_to_calv',
                   '(1 | step_id_)',
                   '(0 + dist_to_cover | id)',
                   '(0 + dist_to_cover:cort_ng_g |id)',
                   '(0 + dist_to_cover:cort_ng_g:d_to_calv |id)')
  
  # Set number of random parameters for mapping
  nvar_parm <- 3
  # Set up model without fitting
  model_form <- suppressWarnings(
    glmmTMB(reformulate(model_covs, response = 'case_'),
            family = poisson(), 
            map = list(theta = factor(c(NA, 1:nvar_parm))),
            data = dat, doFit = F))
  # Set variance of random intercept to 10^6
  model_form$parameters$theta[1] <- log(1e6)
  # Fit model using large fixed variance
  model_fit <- glmmTMB:::fitTMB(model_form) %>%
    broom.mixed::tidy() %>%
    filter(term %in% c('dist_to_cover',
                       'dist_to_cover:cort_ng_g',
                       'dist_to_cover:cort_ng_g:periodpre',
                       'dist_to_cover:cort_ng_g:d_to_calv')) %>%
    select(term, estimate, std.error)

  return(model_fit)
  
}
