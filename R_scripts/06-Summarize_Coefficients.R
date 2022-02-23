
#################################################
###                                           ###
### Calculate bootstrapped model coefficents  ###
### for crop and cover habitats with          ###
### change in cortisol levels (median and     ###
### 95% CI effect sizes)                      ###
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
  for(mods in list.files('alt_output/model_boots/', pattern = i)) {
    # Increase iteration number
    iteration <- iteration + 1
    # Load model
    out_mod <- readRDS(paste('alt_output/model_boots/', mods, sep =''))
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
hab_mod_outs$term <- factor(hab_mod_outs$term, labels = c('cos(TA)', 'Cover', 
                                                          'Crop', 'log(SL)'))
  
# Save model data
saveRDS(all_mods, 'alt_output/model_results_crop-cov.rds')
saveRDS(hab_mod_outs, 'alt_output/model_results_hab.rds')

# Summarize effect sizes (± 95% CI)
mod_summs <- all_mods %>%
  filter(! term %in% c('log(SL)', 'cos(TA)')) %>%
  # Pivot to coefficients in seperate cols
  pivot_wider(names_from = term, values_from = estimate) %>%
  rowwise() %>%
  # Calculate difference in ß across cols for each model
  mutate(general = Habitat,
         post_calving = sum(`Habitat`, `Post:Habitat`),
         general_GC = sum(`Habitat`, `Habitat:GC`),
         post_calving_GC = sum(
           `Habitat` + `Post:Habitat` + `Post:Habitat:GC`)) %>%
  # Pivot back to effect sizes in one col
  pivot_longer(cols = c('general', 'post_calving', 
                             'general_GC', 'post_calving_GC'),
               names_to = 'effect_size') %>%
  # Summarize exp(ß) across models
  group_by(habitat, effect_size) %>%
  summarize(median = exp(median(value)),
            lower = exp(quantile(value, probs = 0.025)),
            upper = exp(quantile(value, probs = 0.975)))
  
# Plot crop-cover models
# tiff('figures/csee_model_bplots_crop-cov.tiff',
#      width = 12, height = 8, units = 'in', res = 300)
ggplot(all_mods, aes(x = term, y = estimate, fill = habitat)) +
  # scale_fill_manual(values = c('#0072e0', '#e06e00')) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_boxplot(width = 0.5, alpha = 0.6, fill = '#b7b4bb') +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        plot.margin = unit(c(0.5, 0.5, 1, 0.5), 'cm'),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 18, colour = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        legend.position = 'none') +
  xlab('') + ylab('Estimate') +
  coord_flip() +
  facet_wrap(~ habitat, scales = 'free_x')
  
  