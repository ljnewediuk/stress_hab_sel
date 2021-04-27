
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

calv_dates <- calv_dates %>%
  rename('animal_year' = animal_ID) %>%
  select(animal_year, calved)

dat <- sample_IDs %>%
  left_join(cort_t3, by = 'label') %>%
  mutate(year = lubridate::year(sample_lmt),
         Jday = lubridate::yday(sample_lmt),
         animal_year = paste(animal_ID, year, sep = '_')) %>%
  left_join(calv_dates, by = 'animal_year') %>%
  mutate(d_to_calv = Jday - calved,
         period = ifelse(d_to_calv < 0, 'pre_calv', 'post_calv'))

horm_means <- dat %>%
  group_by(animal_ID) %>%
  summarize(mean_cort = mean(cort_ng_g, na.rm = T),
            mean_t3 = mean(t3_ng_g, na.rm = T))


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

saveRDS(dat, 'output/horm_calv_dat.rds')


