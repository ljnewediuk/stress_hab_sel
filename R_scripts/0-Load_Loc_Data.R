
#################################################
###                                           ###
### Load and clean location data (if this     ###
### step has not already been completed)      ###
###                                           ###
### Checks input folder and runs cleaning     ###
### if necessary files are not present        ###
###                                           ###
#################################################

# Load packages
library(tidyverse)

# List files in input
list.files('input/')

# If the files are missing, run loop to load and clean them
if(!"vita_elk_lotek_feb_2016-july_2019_cleaned.rds" %in% 
   list.files('input/') &
   !"vita_elk_vectronic_feb_2019-march_2021_cleaned.rds" %in% 
   list.files('input/'))  {
  
  # Source the cleaning function
  source('R_functions/Clean_GPS_Data.R')
  
  # Load data from main files
  system("mkdir ~/Documents/R-Projects/state-dependent_hs/input/temp/")
  system(paste("cp ~/Documents/Elk*Data/Vita*Elk/Collar*data*raw/",
               "vita_elk_lotek_feb_2016-july_2019.csv ~/Documents/R-Projects/",
               "state-dependent_hs/input/temp/", sep=''))
  system(paste("cp ~/Documents/Elk*Data/Vita*Elk/Collar*data*raw/",
               "vita_elk_vectronic_feb_2019-march_2021.csv ~/Documents/R-Projects/",
               "state-dependent_hs/input/temp/", sep=''))
  system(paste("cp ~/Documents/Elk*Data/Vita*Elk/Collar*data*raw/", 
               "collar_deployment_data.csv", 
               " ~/Documents/R-Projects/state-dependent_hs/input/temp/", sep=''))
  
  # Load data into environment
  lotek_dat <- read_csv('input/temp/vita_elk_lotek_feb_2016-july_2019.csv')
  vertex_dat <- read_csv('input/temp/vita_elk_vectronic_feb_2019-march_2021.csv')
  
  # Load collar deployment and retrieval dates into environment
  deploy_dat <- read_csv('input/temp/collar_deployment_data.csv')
  
  # Remove the temp directory
  system("rm -r ~/Documents/R-Projects/state-dependent_hs/input/temp/")
  
  # Join the location data to the deployment dates
  lotek_dat <- lotek_dat %>%
    left_join(deploy_dat, by='animal_ID')
  vertex_dat <- vertex_dat %>%
    left_join(deploy_dat, by='collar_ID')
  
  # Set lotek data to a fix rate of 13 h with 1 h tolerance, maximum elk
  # speed of 45 kph, and WGS84 projection
  lotek_dat <- clean_gps_data(dat = lotek_dat, fix_rate = 13, fix_tol = 1, 
                              fix_unit = 'hours', speed_tol = 45,
                              proj_sys = "+proj=longlat +datum=WGS84 +no_defs")
  
  # Set lotek data to a fix rate of 4 h with 3.8 h tolerance. This is
  # necessary because fix rates vary from 30 min - 4 h; we are just trying to
  # remove bursts of locations ~ 1 sec apart. Also set maximum elk speed of 
  # 45 kph, and WGS84 projection
  vertex_dat <- clean_gps_data(dat = vertex_dat, fix_rate = 4, fix_tol = 3.8, 
                               fix_unit = 'hours', speed_tol = 45,
                               proj_sys = "+proj=longlat +datum=WGS84 +no_defs")
  
  # Additionally remove lotek data from Lac Du Bonnet office (> 50°N)
  lotek_dat <- lotek_dat %>%
    filter(lat < 50)
  
  # Save cleaned data in input folder
  saveRDS(lotek_dat, 'input/vita_elk_lotek_feb_2016-july_2019_cleaned.rds')
  saveRDS(vertex_dat, 'input/vita_elk_vectronic_feb_2019-march_2021_cleaned.rds')
} else {
  cat('Cleaned data already present in input folder')
}

