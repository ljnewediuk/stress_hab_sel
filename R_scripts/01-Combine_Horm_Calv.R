
#################################################
###                                           ###
### Organize hormone and calving data         ###
###                                           ###
### Compile cortisol, T3, calving date data   ###
### in preparation for iSSA models            ###
###                                           ###
### Plot hormones against calving date to     ###
### visualize pattern                         ###
###                                           ###
#################################################

library(tidyverse)

# Load data
sample_IDs <- readRDS('input/final_sample_IDs.rds')
cort_t3 <- read.csv('input/cort_t3_2019-2020.csv')
calv_dates <- readRDS('output/calving_dates.rds')

# Prep calving dates for joining
calv_dates <- calv_dates %>%
  # Rename elk-year column and select only pertinent cols
  rename('animal_year' = animal_ID) %>%
  select(animal_year, calved)

# Join identifications, hormones, and calving dates
dat <- sample_IDs %>%
  left_join(cort_t3, by = 'label') %>%
  mutate(year = lubridate::year(sample_lmt),
         Jday = lubridate::yday(sample_lmt),
         animal_year = paste(animal_ID, year, sep = '_')) %>%
  left_join(calv_dates, by = 'animal_year') %>%
  # Create factor for pre- and post-calving periods
  mutate(d_to_calv = Jday - calved,
         period = ifelse(d_to_calv < 0, 'pre_calv', 'post_calv'))

# Calculate stress response as % increase in over mean cort (baseline)
horm_means <- dat %>%
  group_by(animal_ID) %>%
  mutate(mean_cort = mean(cort_ng_g, na.rm = T),
            mean_t3 = mean(t3_ng_g, na.rm = T)) %>%
  mutate(stress_resp = ifelse(cort_ng_g < mean_cort, 0, cort_ng_g/mean_cort))

# Adjust calv_dates for plotting background
calv_dates <- calv_dates %>%
  # Filter only animals with data
  filter(animal_year %in% unique(horm_means$animal_year)) %>%
  # Add cols expected by ggplot
  mutate(Jday = 150, cort_ng_g = 1)

# Plot cort ad mean cort over calving period by individual
# (remove the two individuals with only 1/2 cort samples (ERE 31 and 20))
# Plot
# tiff('figures/cort_calving_dates.tiff', width = 12, height = 10, units = 'in', res = 300)
ggplot(dat, aes(x = Jday, 
                y = cort_ng_g/1000, 
                group = animal_year)) +
  geom_rect(data = calv_dates,
            aes(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = calved - 3), 
            fill = '#3399ff70') +
  geom_rect(data = calv_dates,
            aes(ymin = -Inf, ymax = Inf, xmin = calved - 3, xmax = calved + 30),
            fill = '#ff7f50', alpha = 0.7) +
  geom_rect(data = calv_dates,
            aes(ymin = -Inf, ymax = Inf, xmin = calved + 30, xmax = Inf),
            fill = '#3399ff70') +
  geom_hline(data = 
               horm_means, 
             aes(yintercept = mean_cort/1000),
             linetype = 'dashed', alpha = 0.6, colour = '#525252') +
  geom_point(size = 3, aes(shape = identification_type)) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', size = 1),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 12, colour = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12, colour = 'black'), 
        axis.title.y = element_text(size = 14, vjust = 5, colour = 'black'),
        axis.title.x = element_text(size = 14, vjust = -5, colour = 'black'),
        legend.position = 'none') +
  ylab('Fecal glucocorticoid metabolites (µ/g)') +
  xlab('Ordinal day') +
  facet_wrap(~ animal_year)

# Save plot
ggsave('figures/MS/cort_calving_dates.tiff', plot = last_plot(), 
       dpi = 300, width = unit(8, 'cm'), height = unit(8, 'cm'), device = 'tiff')

dev.off()

# Save data for remaining prep
saveRDS(horm_means, 'output/horm_calv_dat.rds')


