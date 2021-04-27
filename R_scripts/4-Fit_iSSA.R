
library(sf)
library(tidyverse)
library(amt)

# Make track
trk <- mk_track(all_bursts, 
                .x=X, .y=Y, .t=dat_time, id=animal_ID, cort_ng_g, period,
                crs = sp::CRS("+init=epsg:26914"))

stps <- track_resample(trk, rate=hours(2), tolerance=minutes(30)) %>%
  # Filter to bursts with at least 3 consecutive points (required to calculate TA)
  filter_min_n_burst(min_n=3) %>% 
  steps_by_burst(keep_cols='start') %>% 
  random_steps(n=9) %>%
  group_by(id) %>%
  mutate(log_sl_ = log(sl_ + 1),
         cort_ng_g_sc = scale(cort_ng_g))


# Set max iterations to 10,000 to aid convergence
glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))

# Set up model but don't fit
pop_model <- glmmTMB(case_ ~ log_sl_:cort_ng_g_sc:period + log_sl_ +
                       (1 | step_id_) + 
                       (0 + cort_ng_g_sc | id), # Random effects of step-level density on selection
                     family=poisson(), 
                     data = stps, doFit=FALSE)

# Set variance of random intercept to 10^6
pop_model$parameters$theta[1] <- log(1e6)
nvar_parm <- length(pop_model$parameters$theta)
pop_model$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))

# Fit model using large fixed variance
pop_mod <- glmmTMB:::fitTMB(pop_model)


stps_summ <- stps %>%
  as.data.frame() %>%
  filter(case_ == TRUE) %>%
  group_by(id, cort_ng_g) %>%
  summarize(mean_sl = mean(sl_, na.rm = T))

ggplot(stps_summ, aes(x = cort_ng_g, y = log(mean_sl))) + 
  geom_point() +
  geom_smooth(method = 'lm')




