
#################################################
###                                           ###
### Validate population models using UHC      ###
### plots                                     ###
###                                           ###
#################################################

library(sf)
library(tidyverse)

# Source validation function
source('R_functions/Validate_UHC_RE.R')

# Load data
dat <- readRDS('output/model_dat.rds') %>%
  na.omit() 

# Set seed value
set.seed(1)

# Fit UHC for cover model
cover_uhc <- uhc_validate_re(dat, calv_period = 'pre_calv', 
                             model_form = 'cover', n_iterations = 1000) %>%
  # Filter out covariate for step strata
  filter(! covariate == '(1 | stratum)') %>%
  # Factor covariates
  mutate(covariate = factor(covariate, 
                            labels = c('Cover:FGM', 'Post:Cover:FGM', 'cos(TA)',
                            'Cover', 'Post:Cover', 'log(SL)')))

# UHC for crop model
crop_uhc <- uhc_validate_re(dat, calv_period = 'pre_calv', 
                            model_form = 'crop', n_iterations = 1000) %>%
  # Filter out covariate for step strata
  filter(! covariate == '(1 | stratum)') %>%
  # Factor covariates
  mutate(covariate = factor(covariate, 
                            labels = c('Crop:FGM', 'Post:Crop:FGM', 'cos(TA)',
                                       'Crop', 'Post:Crop', 'log(SL)')))

# Save outputs
saveRDS(cover_uhc, 'output/cover_uhc.rds')
saveRDS(crop_uhc, 'output/crop_uhc.rds')

# Fit and save plots
for(i in c('crop', 'cover')) {
  
  ggplot(get(paste0(i, '_uhc'))) +
    # Plot predicted distribution at used points
    geom_ribbon(aes(x = densdat_x, ymin = densrand_l, ymax = densrand_h), 
                alpha = 0.5) +
    # Plot actual distribution at used points
    geom_line(aes(x = densdat_x, y = densdat_y), colour = 'black', size = 1) +
    # Plot available distribution at used points
    geom_line(aes(x = densavail_x, y = densavail_y), 
              colour = 'red', size = 1, linetype = 'dashed') +
    theme(plot.caption = element_text(size = 22, colour = 'black', hjust = 0.5),
          plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
          axis.title.x = element_text(size = 18, colour = 'black', vjust = -5),
          axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
          panel.background = element_rect(fill = 'white'), 
          panel.grid = element_blank(),
          axis.line.x.bottom = element_line(size = 1, colour = 'black'),
          axis.line.y.left = element_line(size = 1, colour = 'black'),
          axis.text = element_text(size = 16, colour = 'black'),
          strip.text = element_text(size = 16, colour = 'black'),
          strip.placement = 'inside') +
    ylab('Density') +
    xlab('Covariate value') +
    facet_wrap(~ covariate, scales = 'free', strip.position = 'left')
  
  ggsave(paste0('figures/uhc_', i, '.tif'), 
         device = 'tiff', width = 13, height = 10, units = 'in', dpi = 300)
  
  
}

dev.off()
