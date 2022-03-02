
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
  filter(! animal_year %in% 
           c('ER_E_20_2019', 'ER_E_20_2020', 'ER_E_31_2019', 'ER_E_31_2020')) %>%
  # Add cols expected by ggplot
  mutate(Jday = 150, cort_ng_g = 1)

# Plot cort ad mean cort over calving period by individual
# (remove the two individuals with only 1/2 cort samples (ERE 31 and 20))
# Plot
# tiff('figures/cort_calving_dates.tiff', width = 12, height = 10, units = 'in', res = 300)
ggplot(dat[!dat$animal_ID %in% c('ER_E_20', 'ER_E_31') ,], aes(x = Jday, 
                y = cort_ng_g/1000, 
                group = animal_year)) +
  geom_rect(data = calv_dates,
            aes(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = calved), 
            fill = '#5ac18e30') +
  geom_rect(data = calv_dates,
            aes(ymin = -Inf, ymax = Inf, xmin = calved, xmax = Inf),
            fill = '#e4bb3f', alpha = 0.3) +
  geom_point(size = 3, aes(shape = identification_type)) +
  geom_hline(data = 
               horm_means[!horm_means$animal_ID %in% c('ER_E_20', 'ER_E_31') ,], 
             aes(yintercept = mean_cort/1000),
             linetype = 'dashed', alpha = 0.6) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 15),
        panel.grid = element_blank(),
        axis.text = element_text(size = 15), 
        axis.title.y = element_text(size = 18, vjust = 5),
        axis.title.x = element_text(size = 18, vjust = -5),
        legend.position = 'none') +
  ylab('Fecal glucocorticoid metabolites (microgram/g)') +
  xlab('Ordinal day') +
  facet_wrap(~ animal_year)

dev.off()

# Save data for remaining prep
saveRDS(horm_means, 'output/horm_calv_dat.rds')


