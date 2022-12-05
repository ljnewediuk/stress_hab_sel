
#################################################
###                                           ###
### Calculate relative selection strength     ###
### for crop and cover habitats with          ###
### change in cortisol levels, comparing the  ###
### pre- and post-calving periods, using      ###
### bootstrapped model coefficients to get    ###
### median and 95% CI RSS                     ###
###                                           ###
### Then, plot and visualize results          ###
###                                           ###
#################################################

library(tidyverse)

# Load model
mod <- readRDS('output/model_fit.rds')
# Load model data for calculating x ranges
model_dat <- readRDS('derived_data/issa_model_dat.rds') %>%
  # na.omit() %>%
  # Scale  data
  mutate(cort_ng_g_unsc = cort_ng_g,
         dist_to_cover_unsc = dist_to_cover,
         across(c(dist_to_crop, dist_to_cover, cort_ng_g), function(x) as.vector(scale(x))),
         # Make the post-calving period only ~ 30 d when calf most vulnerable
         period = ifelse(d_to_calv > 3 | d_to_calv < -16, 'pre', 'post'),
         # Factor land cover
         across(c(forest, crop, cover), factor),
         # Add log(SL) and cos(TA)
         log_sl_ = log(sl_ + 1),
         cos_ta_ = cos(ta_))

# Set GC quantiles
min_cort <- quantile(model_dat$cort_ng_g, probs = 0.2)
med_cort <- quantile(model_dat$cort_ng_g, probs = 0.5)
max_cort <- quantile(model_dat$cort_ng_g, probs = 0.8)

# Set values for habitat distances
# Crop
med_dist_crop <- median(model_dat$dist_to_crop, na.rm = T)
# Sequence from min to max for cover gradient
cover_quants <- seq(min(model_dat$dist_to_cover), 
                    max(model_dat$dist_to_cover), 
                    length.out = 10)
names(cover_quants) <- paste0(floor(seq(min(model_dat$dist_to_cover_unsc), 
                                        max(model_dat$dist_to_cover_unsc),
                                        length.out = 10)), 'm')
min_dist_cover <- min(model_dat$dist_to_cover, na.rm = T)
max_dist_cover <- max(model_dat$dist_to_cover, na.rm = T)

# Set mean sl_
mean_sl <- mean(model_dat$log_sl_)

# Data frame for predicted data from loc x1
# Maybe instead of max cort, mean of the max corts for all individuals?
x1 <- data.frame(cort_ng_g = seq(from = min_cort, to = max_cort, 
                                 length.out = 100),
                 dist_to_crop = med_dist_crop,
                 step_id_ = NA,
                 log_sl_ = mean_sl,
                 id = NA) 
# Data frame for predicted data from loc x2 (comparison to x1)
x2 <- x1 %>%
  mutate(cort_ng_g = min_cort)

# Loop to calculate RSS

# Initiate dfs to collect results
RSS_preds <- data.frame()
  # Calculate log RSS between loc x1 and x2
  for(per in c('pre', 'post')) {
    for(i in 1:length(cover_quants)) {   
      # Calculate either max or min distance to cover
      # cdist <- ifelse(i == 'max', max(model_dat$dist_to_cover), min(model_dat$dist_to_cover))
      # Add period and distance to cover to data
      x1 <- x1 %>% 
        # mutate(period = per, dist_to_cover = cdist)
        mutate(period = per, dist_to_cover = cover_quants[[i]])
      x2 <- x2 %>% 
        mutate(period = per, dist_to_cover = cover_quants[[i]])
      # Predict and calculate
      logp_1 <- predict(mod, newdata = x1, type = 'link', re.form = NA, se.fit = F)
      logp_2 <- predict(mod, newdata = x2, type = 'link', re.form = NA, se.fit = F)
      # Log RSS
      logRSS <- logp_1 - logp_2
      
      # Assign predicted selection to new df for binding
      RSS_row <- data.frame(id = 1:100,
                            period = per, 
                            log_rss = logRSS,
                            distance_to_cover = names(cover_quants)[[i]], 
                            cort = seq(from = quantile(model_dat$cort_ng_g_unsc, probs = 0.1), 
                                       to = quantile(model_dat$cort_ng_g_unsc, probs = 0.8), 
                                       length.out = 100))
      RSS_preds <- rbind(RSS_preds, RSS_row)

    }
  }

# Pivot into min and max distance cols
RSS_dat <- RSS_preds %>% pivot_wider(names_from = distance_to_cover, values_from = log_rss)

# Save RSS data for plotting
saveRDS(RSS_dat, 'output/rss_preds.rds')


