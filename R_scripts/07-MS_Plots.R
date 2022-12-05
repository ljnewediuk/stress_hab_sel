
#################################################
###                                           ###
### Plot RSS                                  ###
###                                           ###
#################################################

library(tidyverse)
library(glmmTMB)

# Load RSS data
RSS_dat <- readRDS('output/rss_preds.rds')

# Plot RSS gradient from 0 - ~450 m from cover for pre- and post-lactation at 
# different cort levels
ggplot() +
  geom_hline(yintercept = exp(0), linetype = 'dashed', colour = '#cecece') +
  # Add ribbons for gradient
  geom_ribbon(data = RSS_dat, alpha = 0.05,
              aes(x = cort, ymin = exp(`0m`),
                  ymax = exp(`50m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.15,
              aes(x = cort, ymin = exp(`50m`),
                  ymax = exp(`101m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.25,
              aes(x = cort, ymin = exp(`101m`),
                  ymax = exp(`152m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.35,
              aes(x = cort, ymin = exp(`152m`),
                  ymax = exp(`203m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.45,
              aes(x = cort, ymin = exp(`203m`),
                  ymax = exp(`254m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.55,
              aes(x = cort, ymin = exp(`254m`),
                  ymax = exp(`305m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.65,
              aes(x = cort, ymin = exp(`305m`),
                  ymax = exp(`356m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.75,
              aes(x = cort, ymin = exp(`356m`),
                  ymax = exp(`407m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.85,
              aes(x = cort, ymin = exp(`407m`),
                  ymax = exp(`457m`), fill = period)) +
  geom_line(data = RSS_dat, aes(x = cort, y = exp(`457m`), colour = period), size = 1) +
  geom_line(data = RSS_dat, aes(x = cort, y = exp(`0m`), colour = period), linetype = 'dashed') +
  scale_fill_manual(values = c('#ff7f50', '#3399ff')) +
  scale_colour_manual(values = c('#ff4703', '#0073e6')) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        plot.background = element_rect(fill = 'white', colour = NA),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 14, colour = 'black'),
        legend.title = element_text(size = 14, colour = 'black'),
        legend.text = element_text(size = 14, colour = 'black'),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.key = element_rect(fill = 'white'),
        legend.key.size = unit(c(1, 1), 'cm'),
        legend.position = c(0.3, 0.7),
        legend.background = element_rect(fill = 'white'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12, colour = 'black'), 
        axis.ticks = element_line(colour = 'black'),
        axis.title.x = element_text(size = 14, colour = 'black', vjust = -5),
        axis.title.y = element_text(size = 14, colour = 'black', vjust = 5)) +
  xlab('Fecal glucocorticoid metabolite (µg/g)') +
  ylab('log-RSS for distance from cover') 

ggsave('figures/MS/log-RSS.tiff', plot = last_plot(), 
       dpi = 300, width = unit(6, 'cm'), height = unit(6, 'cm'), device = 'tiff')

dev.off()

