
#################################################
### Takes a dataframe of GPS location data    ###
### through several steps of data cleaning:   ###
###                                           ###
###   1. Remove 2D or missing fixes           ###
###   2. Screen out burst fixes               ###
###   3. Keep only locations 1 d after        ###
###      deployment to 1 d before retrieval   ###
###   4. Remove duplicate locations/times     ###
###   5. Remove movements further than        ###
###      possible within fix interval         ###
###                             `             ###`
### Returns a sf object                       ###
###                                           ###
#################################################

# Required columns:
#   animal_ID, collar_ID, time_utc, lat, long, deploy_date, retrieve_date
# Returns columns:
#   animal_ID, collar_ID, lat, long, time_utc, time_lmt, dat_time, geometry 
# Variables:
#   dt: dataframe object, 
#   fix_rate: numeric,
#   fix_unit: 'hours' or 'mins', 
#   fix_tol: lowest acceptable difference from fix rate), 
#   proj_sys ("+proj=longlat +datum=WGS84 +no_defs")
#   speed_tol (kph; 45 for elk)

# Load required packages
library(tidyverse)
library(sf)

# Function
clean_gps_data <- function(dat, fix_rate, fix_tol, fix_unit, proj_sys, speed_tol) {
  
  # Remove row if not a 3D fix
  cat('Removing 2D fixes... \n')
  
  rm_rows <- c()
  for(i in 1:nrow(dat)) {
    if(!grepl('3D', dat[i,]$status)) {
      rm_rows <- c(rm_rows, i)
    }
  }
  
  dat <- dat[ - rm_rows ,]
  
  cat('Rows removed:', length(rm_rows), '\n')
  cat('Done.\n\n')
  
  # Add seconds to datetime to create datetime column
  dat <- dat %>%
    mutate(dat_time = as.POSIXct(paste(dat$time_utc, ':00', sep=''), 
                                 format = "%Y-%m-%d %H:%M:%S", tz='gmt'))
  
  # Screen out burst fixes
  cat('Screening locations to within relocation tolerance... \n')
  
  rm_bursts <- c()
  keep_times <- c()
  for(i in unique(dat$animal_ID)) {
    # Subset by individual, keeping row index and arranging by date
    ID_dat <- dat %>%
      rownames_to_column('index') %>%
      filter(animal_ID == i) %>%
      arrange(dat_time)
    # List fixes outside min tolerance (i.e. bursts)
    for(j in 2:nrow(ID_dat)) {
      time_btwn <- as.numeric(
        difftime(ID_dat[j ,]$dat_time, ID_dat[j - 1 ,]$dat_time, units=fix_unit))
      if(time_btwn < fix_rate-fix_tol) {
        rm_bursts <- c(rm_bursts, as.numeric(ID_dat[j ,]$index))
      } else {
        keep_times <- c(keep_times, time_btwn)
      }
    }
  }
  # Remove fixes outside range
  if(length(rm_bursts) > 0){
    dat <- dat[ - rm_bursts ,]
  }
  
  cat('Rows removed:', length(rm_bursts), '\n')
  cat('Done. ', '(mean relocation rate =', 
      round(mean(keep_times), digits = 0), fix_unit, ';', 'range =',
      round(min(keep_times), digits=2), '-',
      round(max(keep_times), digits=2), fix_unit, 
      ')\n\n', sep = ' ')
  
  # Remove data fom before deployment and after retrieval
  cat('Removing locations before and after deployment... \n')
  
  rm_undeployed <- c()
  for(i in 1:nrow(dat)) {
    if(dat[i,]$dat_time <= as.Date(dat[i,]$deploy_date)+1 | 
       dat[i,]$dat_time >= as.Date(dat[i,]$retrieve_date)-1) {
      rm_undeployed <- c(rm_undeployed, i)
    }
  }
  # Remove locations before and after deployment
  if(length(rm_undeployed) > 0){
    dat <- dat[- rm_undeployed ,]
  }
  
  cat('Rows removed:', length(rm_undeployed), '\n')
  cat('Done.\n\n')
  
  cat('******************************\n', 'Converting to simple features\n\n')
  
  dat <- dat %>% 
    st_as_sf(coords = c('long', 'lat'), remove = F) %>%
    st_set_crs(proj_sys)
  
  # Remove duplicate fixes or timestamps
  cat('Removing duplicate fixes...\n')
  
  dups <- c()
  for(i in unique(dat$animal_ID)) {
    cat('Animal:', i, sep=' ')
    # Subset by individual, keeping row index and arranging by date
    ID_dat <- dat %>%
      rownames_to_column('index') %>%
      filter(animal_ID == i)
    prog_bar <- txtProgressBar(min=2, max=nrow(ID_dat), initial=2)
    for(j in 2:nrow(ID_dat)) {
      setTxtProgressBar(prog_bar, j)
      # Make list of rows with duplicate locations
      if(ID_dat[j ,]$long==ID_dat[j-1 ,]$long & ID_dat[j ,]$lat==ID_dat[j-1 ,]$lat | 
         ID_dat[j ,]$dat_time==ID_dat[j-1 ,]$dat_time) {
        dups <- c(dups, as.numeric(ID_dat[j ,]$index))
      }
      close(prog_bar)
    }
    cat('\n')
  }
  # Remove the duplicates
  if(length(dups) > 0){
    dat <- dat[ - dups ,]
  }
  
  cat('Rows removed:', length(dups), '\n')
  cat('Done.\n\n')
  
  # Filter location points with speed greater than max speed of animal
  cat('Removing speeds greater than ', speed_tol, ' kph... \n', sep='')
  
  fast_locs <- c()
  keep_speeds <- c()
  for(i in unique(dat$animal_ID)) {
    cat('Animal:', i, sep=' ')
    # Subset by individual, keeping row index and arranging by date
    ID_dat <- dat %>%
      rownames_to_column('index') %>%
      filter(animal_ID == i)
    prog_bar <- txtProgressBar(min=2, max=nrow(ID_dat), initial=2)
    # List fixes exceeding speed tolerance
    for(j in 2:nrow(ID_dat)) {
      setTxtProgressBar(prog_bar, j)
      # Calculate interval between fixes
      fix_interval <- abs(as.numeric(difftime(ID_dat[j-1,]$dat_time, 
                                              ID_dat[j,]$dat_time, units='hours')))
      # Calculate geometric distance between consecutive points, 
      # divide by relocation rate, then by 1000 to get kph from metres ph
      if((as.numeric(st_distance(ID_dat[j-1 ,]$geometry, 
                                 ID_dat[j ,]$geometry))/fix_interval)/1000 > speed_tol) {
        fast_locs <- c(fast_locs, as.numeric(ID_dat[j ,]$index))
      } else {
        keep_speeds <- c(keep_speeds, fix_interval/1000)
      }
      close(prog_bar)
    }
    cat('\n')
  }
  # Remove fixes with speed exceeding tolerance
  if(length(fast_locs) > 0){
    dat <- dat[ - fast_locs ,]
  }
  
  
  cat('Rows removed:', length(fast_locs), '\n')
  cat('Done. ', '(mean speed =', 
      round(mean(keep_speeds), digits = 0), 'kph', ';', 'range =',
      round(min(keep_speeds), digits=2), '-',
      round(max(keep_speeds), digits=2), 'kph', 
      ')\n\n', sep = ' ')
  
  # Select columns for final output
  final_dat <- dat %>% select(animal_ID, collar_ID, lat, long, time_utc, time_lmt, dat_time, geometry)
  
  return(final_dat)
  
}
