
library(tidyverse)
library(sf)
library(raster)
library(amt)

# Make directory to house rasters
dir.create('input/temp/')
# Copy into directory
system(paste('cp ~/Documents/Spatial*Data/Manitoba*Data/landcover/ACI/aci_2019.tif',
             '~/Documents/R-Projects/risk_behaviour_repro_copy/input/temp/', sep = ' '))
system(paste('cp ~/Documents/Spatial*Data/Manitoba*Data/landcover/ACI/aci_2020.tif',
             '~/Documents/R-Projects/risk_behaviour_repro_copy/input/temp/', sep = ' '))
# Load rasters
lc_2019 <- raster('input/temp/aci_2019.tif')
lc_2020 <- raster('input/temp/aci_2020.tif')

# Make vectors of crop and cover habitat
forest_hab <- c(210, 220, 230)
crop_hab <- c(133:137, 139, 143, 145:147, 153, 154, 157, 158, 162, 167, 177, 192,
              193, 195, 197)
cover_hab <- c(50, 210, 220, 230)

# Create rasters of crop and cover habitat within 250 m buffer
for (yr in c(2019, 2020)){
  for(i in c('crop', 'cover', 'forest')) {
    
    # Make rasters with either 1/0 binary vals
    bin_rast <- get(paste0('lc_', yr))
    # Make all values of crop/cover = 1, and otherwise zero
    bin_rast[bin_rast %in% get(paste(i, 'hab', sep = '_'))] <- 1
    bin_rast[bin_rast > 1] <- 0
    
    # Make distance-to rasters
    na_rast <- get(paste0('lc_', yr))
    # Make all values either 1 or NA
    na_rast[na_rast %in% get(paste(i, 'hab', sep = '_'))] <- 1
    na_rast[na_rast > 1] <- NA
    dist_rast <- distance(na_rast)
    # Name values (either crop or cover)
    names(dist_rast) <- paste0('dist_to_', i)
    names(bin_rast) <- i
    # Return raster
    assign(paste(i, 'distance_to_raster', yr, sep = '_'), dist_rast)
    assign(paste(i, 'binary_raster', yr, sep = '_'), bin_rast)
    # Save raster
    writeRaster(bin_rast, paste0('rasters/', 
                                 i, 'binary_rast_', yr, '.tif'), overwrite = T)
    writeRaster(dist_rast, paste0('rasters/', 
                                  i, 'distance_to_rast_', yr, '.tif'), overwrite = T)
    
  }
}