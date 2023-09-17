
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
  # Scale distance data
  mutate(across(c(dist_to_cover, cort_ng_g), function(x) as.vector(scale(x))))

# Set seed value
set.seed(1)

# Fit UHC
uhc <- uhc_validate_re(dat, model_spec = 'full', n_iterations = 1000)

# Save outputs
saveRDS(uhc, 'output/uhc_models.rds')

# Fit and save plots
ggplot(uhc[uhc$covariate == 'dist_to_cover' ,]) +
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
    # facet_wrap(~ covariate, scales = 'free', nrow = 2, strip.position = 'left')
  
  ggsave('figures/MS/uhc_plot.tiff', 
         device = 'tiff', width = 8, height = 7, units = 'in', dpi = 300)

dev.off()
