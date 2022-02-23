
#################################################
###                                           ###
### Validate individual models using UHC      ###
### plots                                     ###
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

if(! file.exists('output/uhc_validation_df.rds')) {
  model_val_df <- data.frame()
  
  for(form in c('crop', 'cover')) {
    for(per in c('pre_calv', 'post_calv')) {
      for(i in c(15, 16, 19, 21, 23, 28)) {
        
        try({
          model_val <- uhc_validate(dat = dat, calv_period = per, model_form = form, 
                                    elk_id = paste0('ER_E_', i), n_iterations = 1000)
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
  
  # Save
  saveRDS(uhc_plot_df, 'output/uhc_validation_df.rds')
  
} else {
  
  uhc_plot_df <- readRDS('output/uhc_validation_df.rds')
  
}


uhc_p_list <- list()
combs <- unique(uhc_plot_df[,c('covariate', 'period')])
for(i in 1:nrow(combs)) {
  # Plot
  uhc_p_list[[paste(combs[i ,][[1]], combs[i ,][[2]], sep = '-')]] <- 
    grid::grob(
      ggplot(uhc_plot_df[uhc_plot_df$covariate == combs[i ,][[1]] & 
                           uhc_plot_df$period == combs[i ,][[2]], ]) +
        # Plot predicted distribution at used points
        geom_ribbon(aes(x = densdat_x, ymin = densrand_l, ymax = densrand_h), 
                    alpha = 0.5) +
        # Plot actual distribution at used points
        geom_line(aes(x = densdat_x, y = densdat_y), colour = 'black', size = 1) +
        # Plot available distribution at used points
        geom_line(aes(x = densavail_x, y = densavail_y), 
                  colour = 'red', size = 1, linetype = 'dashed') +
        theme(plot.caption = element_text(size = 22, colour = 'black', hjust = 0.5),
              axis.title = element_blank(),
              panel.background = element_rect(fill = 'white'), 
              panel.grid = element_blank(),
              axis.line.x.bottom = element_line(size = 1, colour = 'black'),
              axis.line.y.left = element_line(size = 1, colour = 'black'),
              axis.text = element_text(size = 16, colour = 'black'),
              strip.text = element_text(size = 16, colour = 'black'),
              strip.placement = 'inside') +
        labs(caption = combs[i ,][[1]]) +
        facet_wrap(~ id, scales = 'free', strip.position = 'left', ncol = 1)
    )
}

y_lab <- grid::textGrob(label = 'Density', rot = 90, gp = grid::gpar(fontsize = 22))
marg <- grid::textGrob(label = '', gp = grid::gpar(fontsize = 20))

# tiff('figures/post-calv_val.tiff', width = 17, height = 10, units = 'in', res = 300)
gridExtra::grid.arrange(y_lab, uhc_p_list$`crop-post_calv`[[1]], 
                        uhc_p_list$`crop:cort_ng_g_sc-post_calv`[[1]],
                        marg, uhc_p_list$`cover-post_calv`[[1]], 
                        uhc_p_list$`cover:cort_ng_g_sc-post_calv`[[1]],
                        ncol = 6, widths = c(0.5, 2, 2, 0.5, 2, 2))

dev.off()

# tiff('figures/pre-calv_val.tiff', width = 17, height = 10, units = 'in', res = 300)
gridExtra::grid.arrange(y_lab, uhc_p_list$`crop-pre_calv`[[1]], 
                        uhc_p_list$`crop:cort_ng_g_sc-pre_calv`[[1]],
                        marg, uhc_p_list$`cover-pre_calv`[[1]], 
                        uhc_p_list$`cover:cort_ng_g_sc-pre_calv`[[1]],
                        ncol = 6, widths = c(0.5, 2, 2, 0.5, 2, 2))

dev.off()
