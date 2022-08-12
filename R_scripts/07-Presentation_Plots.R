
library(tidyverse)
library(glmmTMB)

# Load RSS data
RSS_dat <- readRDS('output/rss_preds.rds')
# Load model
mod <- readRDS('output/model_fit.rds')
# Load hormone data
horm_dat <- read.csv('input/cort_t3_2019-2020.csv')
# Load nmds data
nmds_habs <- readRDS('input/nmds_habs.rds')
nmds_hull <- readRDS('input/nmds_hull.rds')
# Load calving date-Psi data
calv_psi <- readRDS('input/calv_dates_psi.rds')

# Plot boxplots
# Tidy model
mod_coeffs <- broom.mixed::tidy(mod) %>%
  select(term:std.error) %>%
  filter(term %in% c('dist_to_cover', 'dist_to_cover:cort_ng_g')) %>%
  mutate(term = factor(term, labels = c('Distance to cover', 
                                        'Distance to cover:GC')))
# Plot
ggplot() +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = '#e5e5e5') + 
  geom_point(size = 4, colour = '#e5e5e5',
             data = mod_coeffs, aes(x = term, y = estimate)) +
  geom_linerange(size = 1.5, colour = '#e5e5e5',
                 data = mod_coeffs, 
                 aes(x = term, 
                     ymin = estimate - std.error,
                     ymax = estimate + std.error)) +
  theme(panel.background = element_rect(colour = '#e5e5e5', fill = '#36454f'),
        plot.background = element_rect(fill = '#36454f', colour = NA),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 16, colour = '#e5e5e5'), 
        axis.ticks = element_line(colour = '#e5e5e5'),
        axis.title.x = element_text(size = 18, colour = '#e5e5e5', vjust = -5),
        axis.title.y = element_blank()) +
  ylab('Coefficient estimate') 

# Plot RSS gradient from 0 - ~450 m from cover for pre- and post-lactation at 
# different cort levels
# tiff('figures/rss_plot.tiff', width = 10, height = 8, units = 'in', res = 300)
ggplot() +
  geom_hline(yintercept = exp(0), linetype = 'dashed', colour = '#e5e5e5') +
  # Add ribbons for gradient
  geom_ribbon(data = RSS_dat, alpha = 0.1,
              aes(x = cort, ymin = exp(`0m`),
                  ymax = exp(`50m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.2,
              aes(x = cort, ymin = exp(`50m`),
                  ymax = exp(`101m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.3,
              aes(x = cort, ymin = exp(`101m`),
                  ymax = exp(`152m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.4,
              aes(x = cort, ymin = exp(`152m`),
                  ymax = exp(`203m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.5,
              aes(x = cort, ymin = exp(`203m`),
                  ymax = exp(`254m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.6,
              aes(x = cort, ymin = exp(`254m`),
                  ymax = exp(`305m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.7,
              aes(x = cort, ymin = exp(`305m`),
                  ymax = exp(`356m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.8,
              aes(x = cort, ymin = exp(`356m`),
                  ymax = exp(`407m`), fill = period)) +
  geom_ribbon(data = RSS_dat, alpha = 0.9,
              aes(x = cort, ymin = exp(`407m`),
                  ymax = exp(`457m`), fill = period)) +
  scale_fill_manual(values = c('#ff4040', '#3399ff')) +
  theme(panel.background = element_rect(colour = '#e5e5e5', fill = '#36454f'),
        plot.background = element_rect(fill = '#36454f', colour = NA),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        strip.background = element_rect(fill = '#36454f'),
        strip.text = element_text(size = 18, colour = '#e5e5e5'),
        legend.title = element_text(size = 18, colour = '#e5e5e5'),
        legend.text = element_text(size = 18, colour = '#e5e5e5'),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.key = element_rect(fill = '#e5e5e5'),
        legend.key.size = unit(c(1, 1), 'cm'),
        legend.position = c(0.3, 0.7),
        legend.background = element_rect(fill = '#36454f'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 16, colour = '#e5e5e5'), 
        axis.ticks = element_line(colour = '#e5e5e5'),
        axis.title.x = element_text(size = 18, colour = '#e5e5e5', vjust = -5),
        axis.title.y = element_text(size = 18, colour = '#e5e5e5', vjust = 5)) +
  xlab('Fecal glucocorticoid metabolite (microgram/g)') +
  ylab('log-RSS for habitat') 

dev.off()

# Plot GC-T3
ggplot() +
  geom_point(size = 4, colour = '#e5e5e5', alpha = 0.8,
             data = horm_dat, aes(y = cort_ng_g, x = t3_ng_g)) +
  theme(panel.background = element_rect(colour = '#e5e5e5', fill = '#36454f'),
        plot.background = element_rect(fill = '#36454f', colour = NA),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 16, colour = '#e5e5e5'), 
        axis.ticks = element_line(colour = '#e5e5e5'),
        axis.title.x = element_text(size = 18, colour = '#e5e5e5', vjust = -5),
        axis.title.y = element_text(size = 18, colour = '#e5e5e5', vjust = -5))

# Plot NMDS
ggplot() +
  theme(panel.background = element_rect(colour = '#e5e5e5', fill = '#36454f'),
        plot.background = element_rect(fill = '#36454f', colour = NA),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 16, colour = '#e5e5e5'), 
        axis.ticks = element_line(colour = '#e5e5e5'),
        axis.title.x = element_text(size = 18, colour = '#e5e5e5', vjust = -5),
        axis.title.y = element_text(size = 18, colour = '#e5e5e5', vjust = -5)) +
  # Add hulls (by individual)
  scale_colour_viridis_d() +
  geom_polygon(data = nmds_hull, fill = NA,
               aes(x = MDS1, y = MDS2, colour = group)) +
  # Add habitat names
  geom_text(data = nmds_habs, size = 5,
            aes(x = MDS1, y = MDS2), label = rownames(nmds_habs))

# Plot relationship between calving dates and specialization
ggplot(calv_psi, aes(x = calved_unsc, y = psi)) + 
  geom_point() + 
  geom_smooth(method = 'lm', colour = 'black') +
  annotate(geom = 'text', x = 160, y = 1, label = 'R2 = 0.03') +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  xlab('Calving date') + ylab('Specialization index')






