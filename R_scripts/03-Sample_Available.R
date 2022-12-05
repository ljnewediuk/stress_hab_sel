
#################################################
###                                           ###
### Sample available points from MCP and      ###
### steps for iSSA                            ###
###                                           ###
#################################################

library(tidyverse)
library(sf)

# Load data
dat <- readRDS('output/stress_dat.rds') %>% 
  ungroup()

# Make track and resample for iSSA
issa_dat <- dat %>%
  # Add projected coordinates then drop geometry
  mutate(long = st_coordinates(dat)[, 1], lat = st_coordinates(dat)[, 2]) %>%
  st_drop_geometry() %>%
  amt::mk_track(.x = long, .y = lat, .t = dat_time, id = animal_ID, 
                          # Keep data columns
                          uid, cort_ng_g, t3_ng_g, period, d_to_calv,
                          crs = sp::CRS("+init=epsg:26914")) %>%
  # Extract random steps (40 available/used)
  track_resample(rate=minutes(30), tolerance=minutes(5)) %>%
  # Make sure bursts have at least three points
  filter_min_n_burst(min_n = 3) %>% 
  steps_by_burst(keep_cols = 'start') %>%
  random_steps(n = 40)

# Sample available points from MCP for RSFs
rsf_dat <- data.frame()
for(i in unique(dat$uid)) {
  
  # Filter data for individual
  sub_dat <- dat %>%
    filter(uid == i) %>%
    # Drop lat and long
    select(! c(lat, long))
  # Stop if not enough locations
  if(nrow(sub_dat) < 5) next
  # Convert to Spatial Points and get 100% MCP
  sub_sp <- as(sub_dat, 'Spatial')
  sub_mcp <- suppressWarnings(adehabitatHR::mcp(sub_sp, percent = 100))
  # Convert mcp back to sf
  sub_mcp_sf <- as(sub_mcp, 'sf')
  
  # Sample available points (10 x used) and get coords
  sub_avail <- sub_mcp_sf %>%
    st_sample(size = nrow(sub_dat)*10, type = 'random') %>%
    st_coordinates()
  
  # Remove extra columns from used data and bind to available
  sub_used <- sub_dat %>%
    select(! dat_time, jday) %>%
    mutate(case = 1)
  
  # Build df with columns corresponding to used data
  full_dat <- data.frame(long = sub_avail[,1],
                         lat = sub_avail[,2]) %>%
    # Add corresponding columns
    mutate(animal_ID = unique(sub_dat$animal_ID),
           uid = i,
           label = unique(sub_dat$label),
           yr = unique(sub_dat$yr),
           jday = NA,
           d_to_calv = NA,
           calved = NA,
           period = NA,
           cort_ng_g = unique(sub_dat$cort_ng_g),
           t3_ng_g = unique(sub_dat$t3_ng_g),
           case = 0) %>%
    # Convert to sf
    st_as_sf(coords = c('long', 'lat'), crs = st_crs(sub_dat)) %>%
    # Bind to observed data
    rbind(sub_used)
  
  # Bind ua data
  rsf_dat <- rbind(rsf_dat, full_dat)
}

# Save data for extractions
saveRDS(rsf_dat, 'derived_data/rsf_data.rds')
saveRDS(issa_dat, 'derived_data/issa_data.rds')


