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
library(stars)
library(amt)
library(glmmTMB)

# Load data
dat <- readRDS('output/stress_dat.rds') %>%
  ungroup()

# Extract random steps (1x, 5x, 10x, 20x, 30x, 40x 50x, 100x used steps)
avail_samples <- data.frame()

if(file.exists('output/N_avail_df.rds')) {
  
  avail_samples <- readRDS('output/N_avail_df.rds')
  
} else {
  
  for(n_steps in c(1, 5, 10, 20, 30, 40, 50, 100, 500, 1000)) {
    
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
      
      dist_rast <- raster(paste0('rasters/coverdistance_to_rast_', year, '.tif'))
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
    

    
    # Set max iterations to 10^15 to aid convergence
    glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))
    
    # Cover model
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
              data = issa_dat, doFit = F))
    # Set variance of random intercept to 10^6
    model_form$parameters$theta[1] <- log(1e6)
    # Fit model using large fixed variance
    model_fit <- glmmTMB:::fitTMB(model_form) %>%
      broom.mixed::tidy() %>%
      mutate(n_steps) %>%
      filter(term %in% c('dist_to_cover', 'I(log_sl_)',
                         'dist_to_cover:cort_ng_g', 
                         'dist_to_cover:cort_ng_g:periodpre')) %>%
      select(n_steps, term, estimate, std.error)
    
    avail_samples <- rbind(avail_samples, model_fit)
    
  }
  
  # Save
  saveRDS(avail_samples, 'output/N_avail_df.rds')
  
}

# Factor covariates so they appear in order
avail_samples <- avail_samples %>%
  mutate(term_f = factor(term, levels = c('dist_to_cover', 
                       'dist_to_cover:cort_ng_g', 
                       'dist_to_cover:cort_ng_g:periodpre', 'I(log_sl_)'),
                       labels = c('Dist. to cover', 'Dist. to cover:FGM',
                                  'Dist. to cover:FGM:Pre', 'log(SL)')))

# Plot coefficient estimates by used:available points
# tiff('figures/s_select_n_avail.tiff', width = 16, height = 9, units = 'in', res = 300)
ggplot() +
  geom_pointrange(data = avail_samples, aes(x = factor(n_steps), y = estimate,
                                            ymin = estimate - std.error*2,
                                            ymax = estimate + std.error*2,
                                            colour = factor(n_steps))) +
  scale_colour_manual(values = c(rep('black', 5), rep('red', 5))) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 15, colour = 'black'),
        legend.position = 'none',
        panel.grid = element_blank(),
        axis.text = element_text(size = 15, colour = 'black'), 
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -5),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  facet_wrap(~ term_f) +
  xlab('Number of available steps') +
  ylab('Coefficient estimate')

ggsave('figures/MS/s_select_n_avail.tiff', plot = last_plot(), 
       dpi = 300, width = unit(11, 'cm'), height = unit(6, 'cm'), device = 'tiff')

dev.off()

