
#################################################
###                                           ###
### Calculate bootstrapped model coefficents  ###
### for crop and cover habitats with          ###
### change in cortisol levels (median and     ###
### 95% CI RSS)                               ###
###                                           ###
### Then, plot and visualize results          ###
###                                           ###
#################################################

library(tidyverse)
library(ggplot2)

for(i in c('crop', 'cover', 'hab')) {
  
  # Set initial iteration
  iteration <- 0
  # Initiate df
  mod_outs <- data.frame()
  for(mods in list.files('output/model_boots/', pattern = i)) {
    # Increase iteration number
    iteration <- iteration + 1
    # Load model
    out_mod <- readRDS(paste('output/model_boots/', mods, sep =''))
    # Tidy the model
    mod_tidy <- broom.mixed::tidy(out_mod) %>%
      # Remove random effects and intercept term
      filter(effect == 'fixed' & term != '(Intercept)') %>%
      # Add column for iteration and habitat
      mutate(iteration = iteration,
             habitat = i) %>%
      # Remove extra columns
      select(habitat, term, estimate, iteration)
    # Bind together with remaining models
    mod_outs <- rbind(mod_outs, mod_tidy)
  }
  
  # Assign
  assign(paste(i, 'mod_outs', sep = '_'), mod_outs)
  
}

# Bind crop-cover models together
all_mods <- crop_mod_outs %>%
  rbind(cover_mod_outs) %>%
  # Factor coefficients
  mutate(term = factor(term))

# Add factor levels
levels(all_mods$term) <- c('cos(TA)', 'Habitat', 'Habitat:GC', 'Post:Habitat:GC', 
                           'Post:Habitat', 'Habitat', 'Habitat:GC', 
                           'Post:Habitat:GC', 'Post:Habitat', 'log(SL)')
hab_mod_outs$term <- factor(hab_mod_outs$term, labels = c('cos(TA)', 'Cover', 'Crop', 'log(SL)'))
  
# Save model data
saveRDS(all_mods, 'output/model_results_crop-cov.rds')
saveRDS(hab_mod_outs, 'output/model_results_hab.rds')

# Summarize by coefficient (± 95% CI)
model_summs <- all_mods %>%
  group_by(habitat, term) %>%
  summarize(mean_coeff = mean(estimate)) %>%
  mutate(exp_coeff = exp(mean_coeff))

# Plot crop-cover models
# tiff('figures/model_bplots_crop-cov.tiff',
#      width = 12, height = 8, units = 'in', res = 300)
ggplot(all_mods, aes(x = term, y = estimate, fill = habitat)) +
  scale_fill_manual(values = c('#0072e0', '#e06e00')) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_boxplot(width = 0.5, alpha = 0.6) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 18),
        panel.grid = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.position = 'none') +
  xlab('') + ylab('Estimate') +
  coord_flip() +
  facet_wrap(~ habitat, scales = 'free_x')

# Plot habitat model boxplots
# tiff('figures/model_bplots_hab.tiff',
#      width = 8, height = 8, units = 'in', res = 300)
ggplot(hab_mod_outs, aes(x = term, y = estimate, fill = term)) +
  scale_fill_manual(values = c('grey', '#0072e0', '#e06e00', 'grey')) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_boxplot(width = 0.5, alpha = 0.6) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.position = 'none') +
  xlab('') + ylab('Estimate') +
  coord_flip() +
  facet_wrap(~ habitat, scales = 'free_x')



  
  