
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

# Plot cort ad mean cort over calving period by individual
ggplot(dat, aes(x = d_to_calv, 
                y = cort_ng_g, 
                group = animal_ID, 
                col = animal_ID,
       shape = identification_type)) +
  geom_line(size = 1) + 
  geom_point() +
  geom_hline(data = horm_means, 
             aes(col = animal_ID, yintercept = mean_cort),
             linetype = 'dashed', alpha = 0.6) + 
  scale_colour_viridis_d() + 
  facet_wrap(~ animal_ID)

# Save data for remaining prep
saveRDS(horm_means, 'output/horm_calv_dat.rds')


