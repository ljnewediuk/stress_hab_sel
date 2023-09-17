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

# Source required functions
source('R_functions/Test_Steps.R')

# Load data
dat <- readRDS('output/stress_dat.rds') %>%
  ungroup()

# Extract random steps (1x, 5x, 10x, 20x, 30x, 40x 50x, 100x used steps)
avail_samples <- data.frame()
for(i in c(1, 5, 10, 20, 30, 40, 50, 100)) {
  # Get track with desired number of available steps
  steps_dat <- test_steps(dat, n_steps = i) 
  # Fit model to steps
  mod_out <- fit_mod(steps_dat) %>%
    mutate(n_steps = i)
  # Bind with rest of df
  avail_samples <- rbind(avail_samples, mod_out)
  
}

# Factor covariates so they appear in order
avail_samples_f <- avail_samples %>%
  mutate(term_f = factor(term, levels = c('dist_to_cover',
                                          'dist_to_cover:cort_ng_g',
                                          'dist_to_cover:cort_ng_g:periodpre',
                                          'dist_to_cover:cort_ng_g:d_to_calv'),
                       labels = c('Dist. to cover', 'Dist. to cover:FGM',
                                  'Dist. to cover:FGM:Pre', 
                                  'Dist. to cover:FGM:D. to calv.')))

# Plot coefficient estimates by used:available points
ggplot() +
  geom_pointrange(data = avail_samples_f, aes(x = factor(n_steps), y = estimate,
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

