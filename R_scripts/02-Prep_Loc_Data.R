
#################################################
###                                           ###
### Join hormone data with location data      ###
### from previous 18 hours                    ###
###                                           ###
### Huber et al. 2003 JWM found that mean     ###
### time for peak cort after ACTH challenge   ###
### is 18 hours                               ###
###                                           ###
### Data are saved for extracting land cover  ###
###                                           ###
#################################################

library(sf)
library(tidyverse)

# Load loc data
loc_dat <- readRDS('input/vita_elk_vectronic_feb_2019-march_2021_cleaned.rds') %>%
  # Convert lmt to central time and add month
  mutate(time_lmt = lubridate::with_tz(dat_time, tzone = 'Canada/Central'),
         month = lubridate::month(time_lmt)) %>%
  # Filter data from May through August
  filter(month %in% c(5: 8)) %>%
  # Round dates
  mutate(time_lmt = lubridate::round_date(time_lmt, '3 mins'))

# Load sample data
hc_dat <- readRDS('input/final_sample_IDs.rds') %>%
  left_join(read.csv('input/cort_t3_2019-2020.csv')) %>%
  na.omit()

# Load calving dates
calv_dat <- readRDS('derived_data/calving_dates.rds') %>%
  mutate(yr = as.numeric(substr(animal_ID, 9, 12)),
         animal_ID = substr(animal_ID, 1, 7)) %>%
  select(animal_ID, yr, calved)

all_bursts <- data.frame()
# Loop through rows of dat
for(i in 1: nrow(hc_dat)) {
  # Set sample value to match
  sample_time <- hc_dat[i ,]$sample_lmt
  # Get row index of loc dat with local time matching sample time
  # Clear variable
  row_ind <- 0
  row_ind <- which(loc_dat$time_lmt == sample_time & 
                     loc_dat$animal_ID == hc_dat[i ,]$animal_ID)
  # If rows don't match, keep adding half hour until they do, up to 3 hours
  if(length(row_ind) == 0) {
    iter <- 0
    repeat {
      iter <- iter + 1
      sample_time <- sample_time + lubridate::minutes(30)
      row_ind <- which(loc_dat$time_lmt == sample_time &
                         loc_dat$animal_ID == hc_dat[i ,]$animal_ID)
      if(length(row_ind > 0)) break
      if(iter > 6) break
    }
  }
  # Skip to next sample if still no match
  if(is_empty(row_ind)) next
  if(row_ind == 0) next
  # Extract rows of burst up to 18 hours before sample
  row_burst <- loc_dat %>%
    filter(animal_ID == hc_dat[i ,]$animal_ID &
             time_lmt %in% seq(loc_dat[row_ind ,]$time_lmt - lubridate::hours(20), 
                               loc_dat[row_ind ,]$time_lmt, by = '30 min')) %>%
    # Add additional columns
    mutate(yr = lubridate::year(time_lmt),
           jday = lubridate::yday(time_lmt),
           label = hc_dat[i ,]$label,
           cort_ng_g = hc_dat[i ,]$cort_ng_g,
           t3_ng_g = hc_dat[i ,]$t3_ng_g) %>%
    select(! c(collar_ID, time_utc, time_lmt, month)) %>%
    # Reproject CRS to utms
    st_transform(crs = st_crs(26914))
  # Bind together
  all_bursts <- rbind(all_bursts, row_burst)
  
}

# Join location data with sample data and calving dates
hormone_dat <- all_bursts %>% 
  # Join calving dates
  left_join(calv_dat)  %>%
  # Add column for day/night (if location is between sunrise and sunset)
  mutate(period = ifelse(jday < calved, 'pre-calv', 'post-calv')) %>%
  # Add column for pre-/post-calving
  mutate(d_to_calv = calved - jday) %>%
  # Give samples UID by individual sample
  group_by(animal_ID) %>%
  mutate(uid = paste(substr(animal_ID, 1, 7), 
                     yr, cumsum(!duplicated(label)), sep = '-')) 

# Save stress/loc data for extracting land cover
saveRDS(hormone_dat, 'output/stress_dat.rds')

