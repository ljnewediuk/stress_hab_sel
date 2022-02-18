
#################################################
###                                           ###
### Validate models using UHC plots           ###
###                                           ###
### Fit model to only individuals with both   ###
### pre- and post-calving data                ###
###                                           ###
#################################################

library(sf)
library(tidyverse)

# Source validation function
source('R_functions/Validate_UHC.R')

# Load data
dat <- readRDS('output/model_dat.rds') %>%
  na.omit() 

# Set seed value
set.seed(1)

model_val_df <- data.frame()

for(form in c('crop', 'cover')) {
  for(per in c('pre_calv', 'post_calv')) {
    for(i in c(15, 16, 18, 19, 20, 21, 23, 25, 27, 28, 29, 31, 32)) {
      
      try({
        model_val <- uhc_validate(dat = dat, calv_period = per, 
                                  model_form = form, elk_id = paste0('ER_E_', i))
        })
      
      model_val_df <- rbind(model_val_df, model_val)
      
    }
  }  
}


uhc_plot_df <- model_val_df %>%
  na.omit() %>%
  filter(! covariate %in% c('cos_ta_', 'log_sl_')) %>%
  arrange(period) %>%
  mutate(period_id = factor(paste(period, id, sep = ' ')))

# tiff('figures/val_plot.tiff', width = 12, height = 18, units = 'in', res = 300)
uhc_p_list <- list()
for(i in unique(uhc_plot_df$covariate)) {
  
  # Plot
  uhc_p_list[[i]] <- grid::grob(
    ggplot(uhc_plot_df[uhc_plot_df$covariate == i ,]) +
      # Plot predicted distribution at used points
      geom_ribbon(aes(x = densdat_x, ymin = densrand_l, ymax = densrand_h), alpha = 0.5) +
      # Plot actual distribution at used points
      geom_line(aes(x = densdat_x, y = densdat_y), colour = 'black', size = 1) +
      # Plot available distribution at used points
      geom_line(aes(x = densavail_x, y = densavail_y), colour = 'red', size = 1, linetype = 'dashed') +
      facet_wrap(~ period_id, scales = 'free', strip.position = 'left', ncol = 1)
  )

}

grid.arrange(uhc_p_list$crop[[1]], uhc_p_list$`crop:cort_ng_g_sc`[[1]],
            uhc_p_list$cover[[1]], uhc_p_list$`cover:cort_ng_g_sc`[[1]],
            ncol = 4)

dev.off()
