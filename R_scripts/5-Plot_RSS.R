
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

# Load models
crop_model <- readRDS('output/crop_model_results.rds')
cover_model <- readRDS('output/cover_model_results.rds')

# Data frame for predicted data from loc x1
x1 <- data.frame(cort_ng_g_sc = seq(from = 0.4, to = 2.7, length.out = 100),
                      hab = 1,
                      step_id_ = NA,
                      id = NA) 
# Data frame for predicted data from loc x2 (comparison to x1)
x2 <- x1 %>%
  mutate(cort_ng_g_sc = 0.4)

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
    logp_1 <- predict(rss_mod, newdata = x1, type = 'link', re.form = NA)
    logp_2 <- predict(rss_mod, newdata = x2, type = 'link', re.form = NA)
    logRSS <- logp_1 - logp_2
    # Assign predicted selection to new df for binding
    assign(paste(i, j, sep = '_'), 
           data.frame(period = j, 
                      habitat = i, 
                      cort = seq(from = 1.045, to = 6.325, length.out = 100), 
                      selection = logRSS))
  }
  # Bind predicted selection to larger data frame
  all_rss <- all_rss %>%
    rbind(get(paste(i, 'post', sep = '_'))) %>%
    rbind(get(paste(i, 'pre', sep = '_')))
}

# Plot
ggplot(all_rss, aes(x = cort, y = selection, group = period, col = period)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_colour_manual(values = c('#e4bb3f', '#5ac18e')) +
  geom_line(size = 1) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18)) +
  facet_wrap(~ habitat) +
  xlab('Fecal cortisol (microgram/g)') +
  ylab('Selection for crop')

