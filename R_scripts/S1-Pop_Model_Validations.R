
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
dat <- readRDS('derived_data/issa_model_dat.rds') %>%
  # na.omit() %>%
  # Scale distance data
  mutate(across(c(dist_to_cover, cort_ng_g), function(x) as.vector(scale(x))),
         # Make the post-calving period only ~ 30 d when calf most vulnerable
         period = ifelse(d_to_calv > 3 | d_to_calv < -16, 'pre', 'post'),
         # Factor land cover
         across(c(forest, crop, cover), factor),
         # Add log(SL) and cos(TA)
         log_sl_ = log(sl_ + 1),
         cos_ta_ = cos(ta_))

# Set seed value
set.seed(1)

# Fit UHC for each covariate set

uhc_ests <- data.frame()

for(i in c('full', 'cort_only', 'cover_only')) {
  
  uhc <- uhc_validate_re(dat, model_spec = i, n_iterations = 1000) %>%
    # Filter out covariate for step strata
    filter(covariate == 'dist_to_cover')
  
  uhc_ests <- rbind(uhc_ests, uhc)
  
}

# Save outputs
saveRDS(uhc_ests, 'output/uhc_models.rds')

# Fit and save plots
ggplot(uhc_ests[uhc_ests$model_spec == 'full' ,]) +
    # Plot predicted distribution at used points
    geom_ribbon(aes(x = densdat_x, ymin = densrand_l, ymax = densrand_h), 
                alpha = 0.5) +
    # Plot actual distribution at used points
    geom_line(aes(x = densdat_x, y = densdat_y), colour = 'black', size = 0.5) +
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
    xlab('Distance to cover (scaled)')
    # facet_wrap(~ model_spec, scales = 'free', nrow = 1, strip.position = 'left')
  
  ggsave('figures/MS/uhc_plot.tiff', 
         device = 'tiff', width = 6, height = 6, units = 'in', dpi = 300)

dev.off()
