
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
library(ggplot2)

# Load model data for calculating x ranges
model_dat <- readRDS('output/model_dat.rds')

# Set median and max cort from all data
med_cort <- median(model_dat$cort_ng_g_sc, na.rm = T)
max_cort <- max(model_dat$cort_ng_g_sc, na.rm = T)

# Set mean sl_ and cos_ta_
mean_ta <- mean(model_dat$cos_ta_, na.rm = T)
mean_sl <- mean(model_dat$log_sl_)

# Data frame for predicted data from loc x1
# Maybe instead of max cort, mean of the max corts for all individuals?
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
# Separate into cover and crop models
for(i in c('cover', 'crop')) {
  # Initiate df to collect results
  mod_rss <- data.frame()
  # Set initial iteration
  iteration <- 0
  # Loop through each model separately
  for(mods in list.files('output/model_boots/', pattern = i)) {
    # Increase iteration number
    iteration <- iteration + 1
    # Print iteration
    print(iteration)
    # Rename hab column as habitat name
    colnames(x1)[2] <- paste(i)
    colnames(x2)[2] <- paste(i)
    # Get the appropriate model
    rss_mod <- readRDS(paste('output/model_boots/', mods, sep = ''))
    # Calculate log RSS between loc x1 and x2
    for(j in c('pre', 'post')) {
      x1 <- x1 %>% 
        mutate(period = paste(j, 'calv', sep = '_'))
      x2 <- x2 %>% 
        mutate(period = paste(j, 'calv', sep = '_'))
      # Predict and calculate
      logp_1 <- predict(rss_mod, newdata = x1, type = 'link', re.form = NA, se.fit = F)
      logp_2 <- predict(rss_mod, newdata = x2, type = 'link', re.form = NA, se.fit = F)
      # Log RSS
      logRSS <- logp_1 - logp_2
      # Assign predicted selection to new df for binding
      assign(paste(i, j, sep = '_'), 
             data.frame(period = j, 
                        habitat = i, 
                        cort = seq(from = median(model_dat$cort_ng_g, na.rm = T)/1000, 
                                   to = max(model_dat$cort_ng_g, na.rm = T)/1000, 
                                   length.out = 100), 
                        iteration,
                        logRSS))
  
    }
    # Bind predicted selection to larger data frame
    mod_rss <- mod_rss %>% 
      rbind(get(paste(i, 'post', sep = '_'))) %>% 
      rbind(get(paste(i, 'pre', sep = '_')))
    
  }
  
  # Assign
  assign(paste(i, 'mod_rss', sep = '_'), mod_rss)
}

# Bind together
all_rss <- cover_mod_rss %>%
  rbind(crop_mod_rss) %>%
# Calculate median and CI of RSS for bootstrap
  group_by(period, habitat, cort) %>%
  summarize(selection = mean(logRSS),
            lower = quantile(logRSS, probs = 0.05),
            upper = quantile(logRSS, probs = 0.95)) %>%
  filter(! cort > 5.5)

# Save RSS data
saveRDS(all_rss, 'output/rss_results.rds')
saveRDS(crop_mod_rss, 'output/crop_rss.rds')
saveRDS(cover_mod_rss, 'output/cover_rss.rds')

# Plot
# tiff('figures/rss_gc_habitat.tiff', width = 10, height = 8, units = 'in', res = 300)
ggplot(all_rss, aes(x = cort, y = exp(selection), group = period, col = period)) +
  geom_hline(yintercept = exp(0), linetype = 'dashed') +
  scale_colour_manual(values = c('#e4bb3f', '#5ac18e')) +
  scale_fill_manual(values = c('#e4bb3f', '#5ac18e')) +
  geom_ribbon(alpha = 0.3, col = NA,
              aes(ymin = exp(selection) - exp(lower),
                  ymax = exp(selection) + exp(upper), fill = period)) +
  geom_line(size = 1) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        plot.margin = unit(c(0.5, 0, 1, 1), 'cm'),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 18, colour = 'black'),
        legend.title = element_text(size = 18, colour = 'black'),
        legend.text = element_text(size = 18, colour = 'black'),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.key.size = unit(c(1, 1), 'cm'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 16, colour = 'black'), 
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -5),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  facet_wrap(~ habitat, scales = 'free') +
  xlab('Fecal glucocorticoid metabolite (microgram/g)') +
  ylab('log-RSS for habitat')

dev.off()
