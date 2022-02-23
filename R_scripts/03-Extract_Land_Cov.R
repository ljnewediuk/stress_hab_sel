
#################################################
###                                           ###
### Extract random steps, load rasters, and   ###
### extract land cover in preparation for     ###
### iSSA                                      ###
###                                           ###
### Save data for models                      ###
###                                           ###
#################################################

library(tidyverse)
library(sf)
library(raster)
library(amt)

# Load data
dat <- readRDS('output/stress_dat.rds')

# Make directory to house rasters
dir.create('input/temp/')
# Copy into directory
system(paste('cp ~/Documents/Spatial*Data/Manitoba*Data/landcover/ACI/aci_2019.tif',
       '~/Documents/R-Projects/risk_behaviour_repro/input/temp/', sep = ' '))
# Load raster
lc_dat <- raster('input/temp/aci_2019.tif')

# Make vectors of crop and cover habitat
crop_hab <- c(158, 146, 153, 147, 136, 137, 
              167, 157, 154, 145, 197, 133,
              162, 177, 193, 135)
cover_hab <- c(50, 210, 220, 230)

# Create rasters of crop and cover habitat within 250 m buffer
for(i in c('crop', 'cover')) {
  hab_rast <- lc_dat
  # Make all values of crop/cover = 1, and otherwise zero
  hab_rast[hab_rast %in% get(paste(i, 'hab', sep = '_'))] <- 1
  hab_rast[hab_rast > 1] <- 0
  # Buffer values within 250 m buffer
  focal_buffer <- focalWeight(hab_rast, d = 125, type = 'circle')
  buffered_rast <- focal(hab_rast, focal_buffer, na.rm = T, pad = T) 
  # Name values (either crop or cover)
  names(buffered_rast) <- paste(i, 'prop', sep = '_')
  names(hab_rast) <- i
  # Return raster
  assign(paste(i, 'prop_raster', sep = '_'), buffered_rast)
  assign(paste(i, 'binary_raster', sep = '_'), hab_rast)
  # Save binary raster
  raster::writeRaster(hab_rast, paste0('rasters/', i, '_rast.tif'), overwrite = T)
}


# Make track
trk <- amt::mk_track(dat, 
                .x=X, .y=Y, .t=dat_time, id=animal_ID, 
                # Keep cort and stress response data
                cort_ng_g, period, stress_resp,
                crs = sp::CRS("+init=epsg:26914"))

# Extract random steps
stps <- track_resample(trk, rate=hours(2), tolerance=minutes(30)) %>%
  # Make sure bursts have at least three points
  filter_min_n_burst(min_n = 3) %>% 
  steps_by_burst(keep_cols = 'start') %>% 
  random_steps(n = 40) %>%
  # Extract land cover
  extract_covariates(cover_prop_raster, where = 'end') %>%
  extract_covariates(crop_prop_raster, where = 'end') %>%
  extract_covariates(cover_binary_raster, where = 'end') %>%
  extract_covariates(crop_binary_raster, where = 'end') %>%
  group_by(id) %>%
  # Scale cort by ID
  mutate(log_sl_ = log(sl_ + 1),
         cos_ta_ = cos(ta_),
         cort_ng_g_sc = scale(cort_ng_g, scale = T, center = F))

# Remove the temp directory
system("rm -r ~/Documents/R-Projects/risk_behaviour_repro/input/temp/")

# Save data for model
saveRDS(stps, 'output/model_dat.rds')

