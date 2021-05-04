
#################################################
###                                           ###
### Calculate relative selection strength     ###
### for crop and cover habitats with          ###
### change in cortisol levels, comparing the  ###
### pre- and post-calving periods             ###
###                                           ###
### Then, plot and visualize results          ###
###                                           ###
#################################################

library(tidyverse)
library(ggplot2)

# Load models
crop_model <- readRDS('output/crop_model_results.rds')
cover_model <- readRDS('output/cover_model_results.rds')
model_dat <- readRDS('output/model_dat.rds')

# Set median and max cort from all data
med_cort <- median(model_dat$cort_ng_g_sc, na.rm = T)
max_cort <- max(model_dat$cort_ng_g_sc, na.rm = T)

# Set mean sl_ and cos_ta_
mean_ta <- mean(model_dat$cos_ta_, na.rm = T)
mean_sl <- mean(model_dat$log_sl_)

# Data frame for predicted data from loc x1
x1 <- data.frame(cort_ng_g_sc = seq(from = med_cort, to = max_cort, 
                                    length.out = 100),
                      hab = 1,
                      step_id_ = NA,
                      cos_ta_ = mean_ta,
                      log_sl_ = mean_sl,
                      id = NA) 
# Data frame for predicted data from loc x2 (comparison to x1)
x2 <- x1 %>%
  mutate(cort_ng_g_sc = med_cort)

# Loop to calculate RSS
all_rss <- data.frame()
for(i in c('cover', 'crop')) {
  # Rename hab column as habitat name
  colnames(x1)[2] <- paste(i)
  colnames(x2)[2] <- paste(i)
  # Get the appropriate model
  rss_mod <- get(paste(i, 'model', sep = '_'))
  # Calculate log RSS between loc x1 and x2
  for(j in c('pre', 'post')) {
    x1 <- x1 %>% 
             mutate(period = paste(j, 'calv', sep = '_'))
    x2 <- x2 %>% 
             mutate(period = paste(j, 'calv', sep = '_'))
    # Predict and calculate
    logp_1 <- predict(rss_mod, newdata = x1, type = 'link', re.form = NA, se.fit = T)
    logp_2 <- predict(rss_mod, newdata = x2, type = 'link', re.form = NA, se.fit = T)
    # Log RSS
    logRSS <- as.numeric(logp_1$fit) - as.numeric(logp_2$fit)
    # Confidence intervals from model SE
    upperRSS <- (as.numeric(logp_1$fit) + as.numeric(logp_1$se.fit)) - 
      (as.numeric(logp_2$fit) + as.numeric(logp_2$se.fit))
    lowerRSS <- (as.numeric(logp_1$fit) - as.numeric(logp_1$se.fit)) - 
      (as.numeric(logp_2$fit) - as.numeric(logp_2$se.fit))
    # Assign predicted selection to new df for binding
    assign(paste(i, j, sep = '_'), 
           data.frame(period = j, 
                      habitat = i, 
                      cort = seq(from = median(model_dat$cort_ng_g, na.rm = T)/1000, 
                                 to = max(model_dat$cort_ng_g, na.rm = T)/1000, 
                                 length.out = 100), 
                      selection = logRSS,
                      upper = upperRSS,
                      lower = lowerRSS))
  }
  # Bind predicted selection to larger data frame
  all_rss <- all_rss %>%
    rbind(get(paste(i, 'post', sep = '_'))) %>%
    rbind(get(paste(i, 'pre', sep = '_')))
}

# Plot
# tiff('figures/rss_gc__habitat.tiff', width = 10, height = 8, units = 'in', res = 300)
ggplot(all_rss, aes(x = cort, y = exp(selection), group = period, col = period)) +
  geom_hline(yintercept = exp(0), linetype = 'dashed') +
  scale_colour_manual(values = c('#e4bb3f', '#5ac18e')) +
  scale_fill_manual(values = c('#e4bb3f', '#5ac18e')) +
  geom_ribbon(alpha = 0.3, col = NA,
              aes(ymin = exp(selection) - exp(lower), 
                  ymax = exp(selection) + exp(upper), fill = period)) +
  geom_line(size = 1) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18)) +
  facet_wrap(~ habitat) +
  xlab('Fecal cortisol (microgram/g)') +
  ylab('log RSS for habitat')

