
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

# Load data
dat <- readRDS('input/vita_elk_vectronic_feb_2019-march_2021_cleaned.rds')
hc_dat <- readRDS('output/horm_calv_dat.rds')

loc_dat <- dat %>%
  # Convert lmt to central time and add month
  mutate(time_lmt = lubridate::with_tz(dat_time, tzone = 'Canada/Central'),
         month = lubridate::month(time_lmt)) %>%
  # Filter data from May through August
  filter(month %in% c(5: 8)) %>%
  # Round dates
  mutate(time_lmt = lubridate::round_date(time_lmt, '3 mins'))

all_bursts <- data.frame()
# Loop through rows of dat
for(i in 1: nrow(hc_dat)) {
  # Set sample value to match
  sample_time <- hc_dat[i ,]$sample_lmt
  # Get row index of loc dat with local time matching sample time
  row_ind <- which(loc_dat$time_lmt == sample_time & 
                     loc_dat$animal_ID == hc_dat[i ,]$animal_ID)
  # If rows don't match, keep adding half hour until they do
  if(length(row_ind) == 0) {
    repeat {
      sample_time <- sample_time + lubridate::minutes(30)
      row_ind <- which(loc_dat$time_lmt == sample_time & 
                         loc_dat$animal_ID == hc_dat[i ,]$animal_ID)
      if(length(row_ind > 0)) break
    }
  }
  # Extract rows of burst up to 18 hours before sample
  row_burst <- loc_dat %>%
    filter(animal_ID == hc_dat[i ,]$animal_ID &
             time_lmt %in% seq(loc_dat[row_ind ,]$time_lmt - lubridate::hours(18), 
                               loc_dat[row_ind ,]$time_lmt, by = '30 min')) %>%
    # Add columns for hormones and parturition
    mutate(cort_ng_g = hc_dat[i ,]$cort_ng_g,
           stress_resp = hc_dat[i ,]$stress_resp,
           t3_ng_g = hc_dat[i ,]$t3_ng_g,
           d_to_calv = hc_dat[i ,]$d_to_calv,
           period = hc_dat[i ,]$period) %>%
    select(! c(collar_ID:time_lmt, month)) %>%
    # Reproject CRS to utms
    st_transform(crs = st_crs(26914))
  # Bind together
  all_bursts <- rbind(all_bursts, row_burst)
  
}

# Keep X and Y cols
all_bursts <- all_bursts %>%
  cbind(st_coordinates(all_bursts))

# Save stress/loc data for extracting land cover
saveRDS(all_bursts, 'output/stress_dat.rds')

