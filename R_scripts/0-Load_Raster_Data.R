
#################################################
###                                           ###
### Download raster data from Earth Engine    ###
###                                           ###
### This step only needs to be completed      ###
### if the raster data are not already        ###
### downloaded to the system                  ###
###                                           ###
#################################################

library(sf)
library(reticulate)
library(rgee)

# Set up Earth Engine credentials
ee_Initialize()

# Load elk data and bind together
lotek_dat <- readRDS('input/vita_elk_lotek_feb_2016-july_2019_cleaned.rds')
vertex_dat <- readRDS('input/vita_elk_vectronic_feb_2019-march_2021_cleaned.rds')
elk_dat <- rbind(lotek_dat, vertex_dat)

# Transform to projected coordinates
elk_dat_utm <- st_transform(elk_dat, crs = st_crs(26914))

# Create bounding box around the elk data
range_bb <- st_bbox(elk_dat_utm)
# Add 2 km buffer
range_bb <- st_bbox(c(xmin = range_bb[[1]] - 2000, ymin = range_bb[[2]] - 2000, 
                    xmax = range_bb[[3]] + 2000, ymax = range_bb[[4]] + 2000),
                    crs = st_crs(26914))

# Define geometry region for data
geometry <- ee$Geometry$Rectangle(
  coords = as.numeric(range_bb),
  proj = "EPSG:26914",
  geodesic = FALSE
)

# Create temporary output folder to save rasters
system("mkdir ~/Documents/R-Projects/risk_behaviour_repro/output/temp/")

# Import AAFC ACI and USDA CDL rasters from GEE
for(i in c('15', '16', '17', '18', '19')){
  cat('Downloading 20', i, ' CDL and ACI rasters...', '\n\n', sep='')
  # Define image year
  img_aci <- ee$Image(paste('AAFC/ACI/20', i, sep=''))$select('landcover')
  img_cdl <- ee$Image(paste('USDA/NASS/CDL/20', i, sep=''))$select('cropland')
  # Force reprojection
  wgs84 <- ee$Projection('EPSG:26914')
  reproj_aci <- img_aci$reproject(crs=wgs84, scale=30)
  reproj_cdl <- img_cdl$reproject(crs=wgs84, scale=30)
  # Download aci raster
  rast_aci <- ee_as_raster(
    image = reproj_aci,
    region = geometry,
    via = "drive"
  )
  # Download cdl raster
  rast_cdl <- ee_as_raster(
    image = reproj_cdl,
    region = geometry,
    via = "drive"
  )
  # Save to environment
  assign(paste('aci_20', i, sep=''), rast_aci)
  assign(paste('cdl_20', i, sep=''), rast_cdl)
  # Temporarily save to output
  raster::writeRaster(get(paste('aci_20', i, sep='')), paste0('output/temp/aci_20', i, '.tif', sep=''))
  raster::writeRaster(get(paste('cdl_20', i, sep='')), paste0('output/temp/cdl_20', i, '.tif', sep=''))
}

# Copy the temp folder to main files
for(i in list.files(path='output/temp/', pattern='cdl')){
  system(paste("cp ~/Documents/R-Projects/risk_behaviour_repro/output/temp/",
               i, " ~/Documents/Spatial*Data/Minnesota*Data/landcover/CDL/",
               sep=''))
}

for(i in list.files(path='output/temp/', pattern='aci')){
  system(paste("cp ~/Documents/R-Projects/risk_behaviour_repro/output/temp/",
               i, " ~/Documents/Spatial*Data/Manitoba*Data/landcover/ACI/",
               sep=''))
}

# Clean up files on Google Drive
ee_clean_container(name = "rgee_backup", type = "drive")

# Remove the temp directory
system("rm -r ~/Documents/R-Projects/risk_behaviour_repro/output/temp/")

