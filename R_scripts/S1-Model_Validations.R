
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

foo <- uhc_validate(dat = dat, calv_period = 'pre_calv', model_form = 'cover_reduced', elk_id = unique(dat$id))
foo2 <- uhc_validate(dat = dat, calv_period = 'post_calv', model_form = 'cover_reduced', elk_id = unique(dat$id))
foo3 <- foo %>% rbind(foo2)

foo4 <- uhc_validate(dat = dat, calv_period = 'pre_calv', model_form = 'crop', elk_id = unique(dat$id))
foo5 <- uhc_validate(dat = dat, calv_period = 'post_calv', model_form = 'crop', elk_id = unique(dat$id))
foo6 <- foo4 %>% rbind(foo5)

# Plot
ggplot(foo6) +
  # Plot predicted distribution at used points
  geom_ribbon(aes(x = densdat_x, ymin = densrand_l, ymax = densrand_h), alpha = 0.5) +
  # Plot actual distribution at used points
  geom_line(aes(x = densdat_x, y = densdat_y), colour = 'black', size = 1) +
  # Plot available distribution at used points
  geom_line(aes(x = densavail_x, y = densavail_y), colour = 'red', size = 1, linetype = 'dashed') +
  facet_grid(rows = vars(period), cols = vars(covariate), scales = 'free')

# Create UHC plots (uhcdensplot function)
uhcdensplot(densdat = denshats$densdat, 
            densrand = denshats$densrand, 
            includeAvail = TRUE, 
            densavail = denshats$densavail) 
mtext(side=3, line=1,  panlabs1[1], cex=1.2, ad=0)
mtext(outer=F, side=2, line=3, "Density")
mtext(outer=F, side=3, line=1, "Deciduous", cex=1.4)
legend(0.3, 10.5, c("Available", "Used", "Predicted"), 
       lty=c(1, 2, 1), lwd=c(2,3, 10), 
       col=c("Black", "red", "gray"), bty="n", cex=1.8)

