
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
library(stars)
library(amt)

# Load data
issa_dat <- readRDS('derived_data/issa_data.rds')
rsf_dat <- readRDS('derived_data/rsf_data.rds')

# Extract land cover for RSFs
lc_data_rsf <- data.frame()
for(id in unique(rsf_dat$uid)) {
  # Filter individual uid
  sub_dat <- rsf_dat %>%
    filter(uid == id) %>% 
    # Remove duplicate timestamps
    distinct(geometry, .keep_all = T)
  # Extract each land cover type at each point
  for(lc in c('crop', 'cover', 'forest')) {
    # Load raster for appropriate year
    bin_rast <- read_stars(paste0('rasters/', lc, 'binary_rast_', unique(sub_dat$yr), '.tif'))
    dist_rast <- read_stars(paste0('rasters/', lc, 'distance_to_rast_', unique(sub_dat$yr), '.tif'))
    # Extract points
    bin_vals <- st_extract(x = bin_rast, at = sub_dat)
    dist_vals <- st_extract(x = dist_rast, at = sub_dat)
    # Rename to cover type
    colnames(dist_vals)[1] <- paste0('dist_to_', lc)
    colnames(bin_vals)[1] <- lc
    # Join to data
    sub_dat <- st_join(sub_dat, bin_vals) %>% st_join(dist_vals)
  }
  # Bind to df
  lc_data_rsf <- rbind(lc_data_rsf, sub_dat)
}

# Extract land cover for iSSA
lc_data_issa <- data.frame()
for(i in unique(issa_dat$uid)) {
  # Filter individual uid
  sub_dat <- issa_dat %>%
    filter(uid == i)
  # Get year
  yr <- as.numeric(substr(unique(sub_dat$uid), 9, 12))
  
  for(lc in c('crop', 'cover', 'forest')) {
    bin_rast <- raster(paste0('rasters/', lc, 'binary_rast_', yr, '.tif'))
    dist_rast <- raster(paste0('rasters/', lc, 'distance_to_rast_', yr, '.tif'))
    # Rename to cover type
    names(dist_rast)[1] <- paste0('dist_to_', lc)
    names(bin_rast)[1] <- lc
    # Extract land cover
    lc_vals <- sub_dat %>%
      extract_covariates(bin_rast) %>%
      extract_covariates(dist_rast, where = 'end') 
    # Combine data
    sub_dat <- sub_dat %>% left_join(lc_vals)
  }
  lc_data_issa <- rbind(lc_data_issa, sub_dat)
}

# Save data for models
saveRDS(lc_data_issa, 'derived_data/issa_model_dat.rds')
saveRDS(lc_data_rsf, 'derived_data/rsf_model_dat.rds')

