
#################################################
###                                           ###
### Estimate unknown parturition dates using  ###
### known elk calving patterns                ###
###                                           ###
### Uses general method from Marchand et al.  ###
### 2021 with a few important differences:    ###
###                                           ###
### - Uses only the RT100 and distance btwn   ###
###   pts variables (no accelerometer data)   ###
### - Parturition window is 5 days with RT100 ###
###   within a 300 m buffer, because the elk  ###
###   tend to move their calves around within ###
###   a small area over several days          ###
###                                           ###
### Accuracy is ~ 65%                         ###
###                                           ###
#################################################

library(tidyverse)
library(sf)

# Load elk data
dat <- readRDS('vita_elk_vectronic_feb_2019-march_2021_cleaned.rds')

# Reproject into UTMs and add as column to data frame
utm_conv <- dat %>% st_transform(crs = 26914) 

dat_utm <- unlist(st_geometry(utm_conv)) %>% 
  matrix(ncol=2,byrow=TRUE) %>% 
  as.data.frame() %>%
  setNames(c('x', 'y')) %>%
  cbind(utm_conv)

# Subset data to calving period and add year to ID
dat_calving <- dat_utm %>%
  mutate(month = lubridate::month(dat_time)) %>%
  filter(month %in% c(5:7)) %>%
  mutate(year = lubridate::year(dat_time + lubridate::hours(6))) %>%
  mutate(animal_ID = paste(animal_ID, year, sep= '_'))

# Format data for recurse
dat_recurse <- dat_calving %>%
  rename(id = animal_ID, t = dat_time) %>%
  mutate(id = as.factor(id)) %>%
  select(id, x, y, t)

# Make df of estimated calving dates
calv_dates <- data.frame( 
  id = c('ER_E_32_2019', 'ER_E_18_2019', 'ER_E_23_2019', 'ER_E_29_2019',
         'ER_E_25_2019', 'ER_E_15_2019', 'ER_E_28_2019', 'ER_E_27_2019',
         'ER_E_24_2019', 'ER_E_30_2019', 'ER_E_26_2019', 'ER_E_20_2019', 
         'ER_E_16_2020', 'ER_E_25_2020', 'ER_E_19_2020', 'ER_E_15_2020'),
  calved = c(135, 141, 143, 144, 151, 153, 154, 159, 159, 161, 162, 163,
             160, 144, 174, 148))
  # Adjusted dates based on graphs
  # calved = c(137, 141, 145, 144, 149, 155, 156, 159, 159, 161, 164, 161,
  #            162, 142, 174, 148))

### Setting up data:

# 1 - recurse::getRecursions for each location point in df (by animal_ID, gives
#     number of locations within specified radius (100 m) of location point;
#     maybe set threshold = 168 h to exclude return visits to the radius beyond
#     one week of potential parturition date)

# Interesting that there seems to be  a pattern, but this seems to be fewest 
# revisits at time of calving. Maybe this is because the elk is not actually 
# leaving the 100 m radius? This is puzzling but maybe the lack of returns will
# be sufficient to pick up the pattern?

# Count recursive visits within 300 m of each point (RT100); wider than in 
# Marchand et al. because elk move the calf between different rest areas
revis <- data.frame()
for(i in unique(dat_recurse$id)){
  
  sub_recurse <- dat_recurse %>%
    filter(id == i)
  
  recursions <- recurse::getRecursions(
  sub_recurse, radius = 300)
  
  id_revis <- data.frame(id = i, t = sub_recurse$t, RT100 = recursions$revisits)
  revis <- rbind(revis, id_revis)
}


# Add revisits and calving dates
dat_recurse_calving <- dat_recurse %>%
  # Add revisit data
  left_join(revis, by = c('id', 't')) %>%
  # Convert t to local time so calving dates are correct and add day/year col
  mutate(time_lmt = t + lubridate::hours(6)) %>%
  mutate(day = lubridate::yday(time_lmt)) 

# Join calving date data in new df
recurse_by_day <- dat_recurse_calving %>%
  select(id, RT100, time_lmt, day) %>%
  right_join(calv_dates, by = 'id') %>%
  # Add column for parturition period (72 h)
  mutate(days_to_calv = day - calved) %>%
  filter(days_to_calv >= -30 & days_to_calv <= 30) %>%
  mutate(parturition = factor(ifelse(days_to_calv %in% seq(0, 5, 1), 'yes', 'no')))

# 2 - Calculate distance between all t and t-1 points for each individual

# Add year to ID column
dat_distance <- dat_utm %>%
  mutate(id = paste(animal_ID, lubridate::year(dat_time), sep = '_')) %>%
  mutate(time_lmt = dat_time + lubridate::hours(6)) %>%
  select(id, x, y, time_lmt)

# Initialize data frame
distance_by_day <- data.frame()

for(i in unique(dat_distance$id)) {
  # Subset by individual id and add distance column
  id_sub <- dat_distance %>%
    filter(id %in% i) %>%
    mutate(dist_last_pt = NA)
  # Loop through rows to calculate distance to previous point
  for(j in 1:nrow(id_sub)) {
    # First row == NA
    if(j == 1) {
      id_sub[j ,]$dist_last_pt <- NA
      # Otherwise, calculate distance between current and previous point
    } else {
      pt_t <- c(id_sub[j ,]$x, id_sub[j ,]$y)
      pt_t_last <- c(id_sub[j-1 ,]$x, id_sub[j-1 ,]$y)
      id_sub[j ,]$dist_last_pt <- raster::pointDistance(pt_t, pt_t_last, lonlat = F)
    }
  }
  # Bind together
  distance_by_day <- rbind(distance_by_day, id_sub)
}

# Join data for model
rf_dat <- left_join(recurse_by_day, distance_by_day, by = c('id', 'time_lmt')) %>%
  na.omit()

# 3 - Plot to visualize data

# Plot RT100 by days to parturition
ggplot(rf_dat, aes(x = days_to_calv, y = RT100, group = id, col = id)) +
  scale_colour_viridis_d() +
  scale_x_continuous(breaks = seq(-30, 30, 5)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 5, linetype = 'dashed') +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf), 
            alpha = 0.2, colour = NA, fill = 'grey') +
  geom_line() + 
  geom_smooth(method = 'loess') + 
  facet_wrap(~id)

# Plot distance between points by days to parturition
ggplot(rf_dat, aes(x = days_to_calv, y = dist_last_pt, group = id, col = id)) +
  scale_colour_viridis_d() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 6, linetype = 'dashed') +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf), 
            alpha = 0.2, colour = NA, fill = 'grey') +
  geom_line() + 
  geom_smooth(method = 'loess') + 
  facet_wrap(~id)

### Fitting model:

# 4 - Run RF model for training data (parturition date ~ RT100 + distance)

# Perform down-sampling
dsample_recipe <- recipes::recipe(parturition~., data = rf_dat) %>% 
  themis::step_downsample(parturition, under_ratio = 1, id = id)
# Prep
dsample_prep <- recipes::prep(dsample_recipe)
# Juice (return pre-processed data)
dsample_rf_dat <- as.data.frame(recipes::juice(dsample_prep))
# Model (without distance variable for now)
rf_mod <- randomForest::randomForest(parturition ~ RT100, data = dsample_rf_dat)
res <- as.data.frame(rfPermute::confusionMatrix(rf_mod, threshold=0.5))

# 5 - Calculate probability of each training data point falling within the 
#     parturition period using predict, type = probability; this returns
#     a two-column matrix with the probability of each location point belonging
#     to either time period (parturition/not)

rf_train_results <- rf_dat %>%
  mutate(part_prob = predict(rf_mod, newdata = rf_dat, type = 'prob')[,2])

### Estimating parturition dates:

# 6 - Set threshold as 1% quantile of parturition = T probabilities in training 
#     data

thresh <- quantile( 
  rf_dat[rf_dat$parturition == 'yes' ,]$part_prob, probs = 0.1, na.rm = T)

# 6 - Calculate probability of each testing data point falling within the 
#     parturition period using the RF model

# Join data frame for elk not included in training data
rf_test_dat <- dat_recurse_calving %>%
  select(id, RT100, time_lmt) %>%
  left_join(distance_by_day, by = c('id', 'time_lmt')) %>%
  # Remove elk already in training set
  filter(!id %in% unique(rf_dat$id))

# 7 - Use zoo::rollapply to calculate mean probability of parturition within 
#     week-long moving window in training set (known elk probability of calving 
#     remains high for approximately 5 days)

# Predict probability of parturition
rf_results <- rf_test_dat %>%
  mutate(part_prob = 
           predict(rf_mod, newdata = rf_test_dat, type = 'prob')[,2]) %>%
  mutate(part_prob_roll = 
           zoo::rollapply(part_prob, width = 120, mean, fill = T)) %>%
  # filter probabilities greater than threshold
  filter(part_prob_roll != 1) %>%
  # Add day col
  mutate(day = lubridate::yday(time_lmt)) %>%
  # factor ID column
  mutate(id = factor(id))

# 8 - If rolling window mean > threshold, identify as probable parturient and 
#     define predicted parturition date/time as peak of smoothed curve;
#     If rolling window mean < threshold for an individual's entire data, 
#     identify that individual as non-parturient

# Visualize the parturition dates for the known elk
ggplot(rf_train_results, 
       aes(x = days_to_calv, y = part_prob, group = id, col = id)) +
  scale_colour_viridis_d() +
  geom_hline(yintercept = thresh, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 5, linetype = 'dashed') +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf), 
            alpha = 0.2, colour = NA, fill = 'grey') +
  geom_smooth(method = 'loess') + 
  facet_wrap(~id)

# Plot the probability of parturition for training elk
rf_test_plot <- ggplot(rf_results, 
                       aes(x = day, y = part_prob_roll, group = id, col = id)) +
  scale_colour_viridis_d() +
  geom_hline(yintercept = thresh, linetype = 'dashed') +
  geom_line() + 
  geom_smooth(method = 'loess', span = 0.5) + 
  facet_wrap(~id)

# Create data frame of ids by plot panel for joining
panel_IDs <- data.frame(PANEL = factor(seq(1, 18, 1)), id = levels(rf_results$id))

# Create a ggplot build object to get maximum of smoothed line
rf_test_build <- ggplot_build(rf_test_plot)

# Make table of possible calving dates by ID
x_at_max_y <- ggplot_build(rf_test_plot)$data[[3]] %>%
  left_join(panel_IDs, by = 'PANEL') %>%
  group_by(PANEL) %>%
  # Subset Julian days between May 15 and July 20 to get rid of edge effects
  filter(x > 135 & x < 201) %>%
  filter(y == max(y))

# Visualize the possible parturition dates for the unknown elk
rf_test_plot + 
  geom_vline(data = x_at_max_y, aes(xintercept = x))

# 9 - Join predicted parturition dates with known parturition dates and save

# Join together known and predicted calving dates
final_calv_dates <- x_at_max_y %>%
  ungroup() %>%
  mutate(x = round(x, digits = 0), predicted = 'yes') %>%
  rename('calved' = x) %>%
  select(id, calved, predicted) %>%
  plyr::rbind.fill(calv_dates) %>%
  mutate(predicted = replace_na(predicted, 'no')) %>%
  rename('animal_ID' = id)

# Save
saveRDS(final_calv_dates, 'output/calving_dates.rds')

