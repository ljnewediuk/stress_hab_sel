
#################################################
###                                           ###
### Calculate relative selection strength     ###
### for distance to cover with                ###
### change in cortisol levels, comparing the  ###
### pre- and post-calving periods and days    ###
### since calving                             ###      
###                                           ###
#################################################

library(tidyverse)

# Load models
mod <- readRDS('output/model_fit_bs.rds')

# Load model data for calculating x ranges
model_dat <- readRDS('derived_data/issa_model_dat.rds') %>%
  # na.omit() %>%
  # Scale  data
  mutate(cort_ng_g_unsc = cort_ng_g,
         dist_to_cover_unsc = dist_to_cover,
         across(c(dist_to_crop, dist_to_cover, cort_ng_g), function(x) as.vector(scale(x))),
         # Add log(SL) and cos(TA)
         log_sl_ = log(sl_ + 1),
         cos_ta_ = cos(ta_))

# Set GC quantiles
min_cort <- quantile(model_dat$cort_ng_g, probs = 0.2)
med_cort <- quantile(model_dat$cort_ng_g, probs = 0.5)
max_cort <- quantile(model_dat$cort_ng_g, probs = 0.8)

# Set values for habitat distances
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
x1 <- data.frame(cort_ng_g = seq(from = min_cort, to = max_cort, 
                                 length.out = 100),
                 step_id_ = NA,
                 log_sl_ = mean_sl,
                 id = NA) 
# Data frame for predicted data from loc x2 (comparison to x1)
x2 <- x1 %>%
  mutate(cort_ng_g = min_cort)

# Function to calculate RSS for all distances from cover
calc_RSS <- function(n_days, per, mod) {
  
  RSS <- data.frame()
  for(i in 1:length(cover_quants)) {   
    # Add period, days to calving, and distance to cover to data
    x1 <- x1 %>% 
      mutate(d_to_calv = n_days, period = per, dist_to_cover = cover_quants[[i]])
    x2 <- x2 %>% 
      mutate(d_to_calv = n_days, period = per, dist_to_cover = cover_quants[[i]])
    # Predict and calculate
    logp_1 <- predict(mod, newdata = x1, type = 'link', re.form = NA, se.fit = F)
    logp_2 <- predict(mod, newdata = x2, type = 'link', re.form = NA, se.fit = F)
    # Log RSS
    logRSS <- logp_1 - logp_2
    
    # Assign predicted selection to new df for binding
    RSS_row <- data.frame(id = 1:100,
                      d_to_calv = n_days, 
                      period = per,
                      log_rss = logRSS,
                      distance_to_cover = names(cover_quants)[[i]], 
                      cort = seq(from = quantile(model_dat$cort_ng_g_unsc, probs = 0.1), 
                                 to = quantile(model_dat$cort_ng_g_unsc, probs = 0.8), 
                                 length.out = 100))
    # Bind together
    RSS <- rbind(RSS, RSS_row)
  }
  
  return(RSS)
  
}
  
# Predict RSS for days since calving, setting period to post-calving
preds_d <- data.frame()
for(x in 1:length(mod)) {
  for(i in c(0, 15, 30, 45, 60)) {
    
    preds_RSS <- calc_RSS(i, 'post-calv', mod[[x]]) %>%
      mutate(mod_iteration = x)
    
    preds_d <- rbind(preds_d, preds_RSS)
  }
}

# Predict RSS for all models pre- and post-calving, setting days since calving to zero
preds_p <- data.frame()
for(x in 1:length(mod)) {
  for(i in c('pre-calv', 'post-calv')) {
    
    preds_RSS <- calc_RSS(0, i, mod[[x]]) %>%
      mutate(mod_iteration = x)
    
    preds_p <- rbind(preds_p, preds_RSS)
    
  }
}

# Save RSS data for plotting (period and days since)
saveRDS(preds_p, 'output/rss_preds_p.rds')
saveRDS(preds_d, 'output/rss_preds_d.rds')

