
#################################################
###                                           ###
### Plots for MS                              ###
###                                           ###
#################################################

library(tidyverse)
library(tidybayes)
library(emmeans)
library(modelr)
library(rstanarm)
library(cowplot)

# Load RSS data
# Period
RSS_p <- readRDS('output/rss_preds_p.rds') %>%
  group_by(period, distance_to_cover, cort) %>%
  summarize(mean_rss = mean(log_rss),
            ci_rss = (sd(log_rss)/sqrt(n()))*2) %>%
  filter(distance_to_cover == '203m')
# Days since calving
RSS_d <- readRDS('output/rss_preds_d.rds') %>%
  group_by(period, d_to_calv, distance_to_cover, cort) %>%
  summarize(mean_rss = mean(log_rss),
            ci_rss = (sd(log_rss)/sqrt(n()))*2) %>%
  filter(distance_to_cover == '203m') %>%
  filter(d_to_calv %in% c(0, 30, 60))

# Load draws from brms model
brms_draws <- readRDS('output/cort_model_draws.rds') %>%
  mutate(period = factor(period, levels = c('post_calv', 'pre_calv'), labels = c('Post', 'Pre')))

# Load hormone/calving dates data
horm_calv_dat <- readRDS('output/horm_calv_dat.rds') %>%
  mutate(period = factor(period, levels = c('post_calv', 'pre_calv'), labels = c('Post', 'Pre')))

# Plot RSS for 200 m from cover pre- and post-lactation at different cort levels
period_RSS <- ggplot() +
  geom_hline(yintercept = exp(0), linetype = 'dashed', colour = '#cecece') +
  geom_smooth(data = RSS_p, se = F, aes(x = cort, y = exp(mean_rss), linetype = period), colour = '#808080') +
  
  geom_ribbon(data = RSS_p, aes(x = cort, ymin = exp(mean_rss) - exp(ci_rss), ymax = exp(mean_rss) + exp(ci_rss), fill = period), alpha = 0.3) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', size = 1),
        plot.background = element_rect(fill = 'white', colour = NA),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        legend.title = element_text(size = 16, colour = 'black'),
        legend.text = element_text(size = 14, colour = 'black'),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.key = element_rect(fill = 'white'),
        legend.key.size = unit(c(1, 1), 'cm'),
        legend.position = c(0.35, 0.85),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 14, colour = 'black'), 
        axis.ticks = element_line(colour = 'black'),
        axis.title.x = element_text(size = 16, colour = 'black', vjust = -5),
        axis.title.y = element_text(size = 16, colour = 'black', vjust = 5)) +
  scale_linetype_manual(values = c('dashed', 'solid')) +
  scale_fill_manual(values = c('#808080', '#808080')) +
  scale_y_continuous(breaks = c(0, 2.5, 5, 7.5)) +
  ylim(-0.2, 9) +
  xlab('Fecal glucocorticoid metabolite (µg/g)') +
  ylab('log-RSS for 200 m from forest/shrub') 

# Plot RSS for 200m from cover at different cort levels and days since calving
days_RSS <- ggplot() +
  geom_hline(yintercept = exp(0), linetype = 'dashed', colour = '#cecece') +
  geom_smooth(data = RSS_d, se = F, aes(x = cort, y = exp(mean_rss), colour = factor(d_to_calv))) +
  
  geom_ribbon(data = RSS_d, aes(x = cort, ymin = exp(mean_rss) - exp(ci_rss), ymax = exp(mean_rss) + exp(ci_rss), fill = factor(d_to_calv)), alpha = 0.3) +
  scale_color_manual(values = c('#ecd0b8', '#ae6134', '#4b2c19'), name = 'Calf age (days)') +
  scale_fill_manual(values = c('#ecd0b8', '#ae6134', '#4b2c19'), name = 'Calf age (days)') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', size = 1),
        plot.background = element_rect(fill = 'white', colour = NA),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        legend.title = element_text(size = 16, colour = 'black'),
        legend.text = element_text(size = 14, colour = 'black'),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.key = element_rect(fill = 'white'),
        legend.key.size = unit(c(1, 1), 'cm'),
        legend.position = c(0.22, 0.82),
        legend.background = element_rect(fill = 'white'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 14, colour = 'black'), 
        axis.ticks = element_line(colour = 'black'),
        axis.title.x = element_text(size = 16, colour = 'black', vjust = -5),
        axis.title.y = element_text(size = 16, colour = 'white', vjust = 5)) +
  scale_y_continuous(breaks = c(0, 2.5, 5, 7.5)) +
  ylim(-0.2, 9) +
  xlab('Fecal glucocorticoid metabolite (µg/g)') +
  ylab('log-RSS for 200 m from forest/shrub') 

# Arrange panels with cowplot
plot_grid(period_RSS, days_RSS, labels = c('A', 'B'), label_size = 18, align = 'hv')

# Save plot
ggsave('figures/MS/log-RSS.tiff', plot = last_plot(), 
       dpi = 300, width = unit(12, 'cm'), height = unit(6, 'cm'), device = 'tiff')

# Plot brms model testing difference in glucocorticoid production pre- and post- calving

ggplot(brms_draws, aes(x = estimate_unsc, y = period)) +
  stat_halfeye(colour = 'black', fill = '#808080', slab_color = 'black', slab_size = 0.5, alpha = 1, slab_alpha = 0.7, linewidth = 5, linetype = 'solid', size = 10, orientation = 'horizontal') +
  geom_jitter(data = horm_calv_dat, aes(x = cort_ng_g, y = period), colour = '#808080', fill = '#808080', pch = 21, alpha = 0.7) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', size = 1),
        plot.background = element_rect(fill = 'white', colour = NA),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        legend.position = 'none',
        panel.grid = element_blank(),
        axis.text = element_text(size = 14, colour = 'black'), 
        axis.ticks = element_line(colour = 'black'),
        axis.title.x = element_text(size = 16, colour = 'black', vjust = -5),
        axis.title.y = element_text(size = 16, colour = 'black', vjust = 5)) +
  ylab('Period') + xlab('Fecal glucocorticoid metabolite (µg/g)')

# Save plot
ggsave('figures/MS/brms_plot.tiff', plot = last_plot(), 
       dpi = 300, width = unit(6, 'cm'), height = unit(6, 'cm'), device = 'tiff')
