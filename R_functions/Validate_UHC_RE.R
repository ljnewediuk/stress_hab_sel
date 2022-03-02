
#################################################
###                                           ###
### Validates model by comparing              ###
### distributions of explanatory variables at ###
### observed and predicted presence locations ###
### (habitat characteristics associated with  ###
### used locations)                           ###
###                                           ###
### Returns df of used/available              ### 
### distributions and 95% CI bounds of        ### 
### predicted habitat characteristics of used ###
### points.                                   ###
###                                           ###
### Function for random effects models fit    ###
### using glmmTMB.                            ###
###                                           ###
### Uses UHC plot process by Fieberg et al.   ###
### 2018  https://doi.org/10.1111/ecog.03123  ###
###                                           ###
#################################################

library(sf)
library(tidyverse)
library(uhcplots)

# Required columns:
#   id, period, case_, step_id_, cover, crop, cort_ng_g_sc, log_sl_, cos_ta_
# Arguments:
#   - Data: Name of data frame
#   - Period: "pre_calv" or "post_calv"
#   - Formula: 
#       - "crop": Contains all crop covariates and crop:cort interaction
#       - "cover": Contains all cover covariates and cover:cort interaction
#       - "crop_reduced": Crop model without interaction
#       - "cover_reduced": Cover model without interaction

uhc_validate_re <- function(dat, calv_period, model_form, n_iterations) {
  
  # Name covariates
  if(model_form == 'cover') covs <- c('cort_ng_g_sc:cover:period', 'cort_ng_g_sc:cover',
                                      'cover:period', 'cover', 'log_sl_', 'cos_ta_',
                                      '(1 | stratum)', '(0 + cort_ng_g_sc + cover | id)')
  if(model_form == 'crop') covs <- c('cort_ng_g_sc:crop:period', 'cort_ng_g_sc:crop', 
                                     'crop:period', 'crop', 'log_sl_', 'cos_ta_',
                                     '(1 | stratum)', '(0 + cort_ng_g_sc + crop | id)')
  # Name of covariates including main effects only
  covs_no_ranef <- covs[1:(length(covs)-2)]
  
  # Rename cases and strata for function variables, add id-stratum for sampling
  indiv_dat <- dat %>%
    rename('stratum' = step_id_) %>%
    mutate(presence = ifelse(case_ == TRUE, 1, 0))
  
  # Separate into training/testing data stratified by individuals/strata
  train.steps <- sample(unique(indiv_dat$stratum), ceiling(length(unique(indiv_dat$stratum))/2))
  test.steps <- unique(indiv_dat$stratum)[! unique(indiv_dat$stratum) %in% train.steps]
  mdat.train <- indiv_dat %>%
    filter(stratum %in% train.steps)
  mdat.test <- indiv_dat %>%
    filter(stratum %in% test.steps)
  
  # Specify model formula
  form1a <- reformulate(covs, response = 'presence')
  # Specify covariates for sampling from multivariate normal distribution
  form2a <- reformulate(c(covs_no_ranef, -1))
  
  # Fit training model
  # Attach glmmTMB
  library(glmmTMB)
  # Set N iterations
  glmmTMBControl(optCtrl=list(iter.max=1e15,eval.max=1e15))
  # Set up model without fitting
  ssf.train <- suppressWarnings(
    glmmTMB(form1a,
            family = poisson(), 
            map = list(theta = factor(c(NA, 1:3))),
            data = mdat.train, doFit = F))
  # Set variance of random intercept to 10^6
  ssf.train$parameters$theta[1] <- log(1e6)
  # Fit model using large fixed variance
  ssf.train <- glmmTMB:::fitTMB(ssf.train)
  
  # Design matrix from test data for SSF
  design.mat.test <- model.matrix(form2a, data = mdat.test)
  # Design matrix for covariates z (matrix of used & available environmental
  # characteristics in the test data; can differ/be the same as design.mat.test)
  z <- model.matrix(form2a, data = mdat.test)
  
  # Create simulation envelopes for the environmental characteristics at the 
  # observed locations in the test data (simstrat function)

  # Define simstrat function for random effects model
  simstrat <- function(nsims, xmat, stratum, fit_ssf, z){
    ustrat <- unique(stratum) # unique strata
    nstrat <- length(ustrat) # number of strata
    # array to store chosen x's for each simulation
    x_sim_choice <- array(NA, dim = c(nsims, nstrat, ncol(z)))
    # new beta^ for each simulation
    # Get coefficients from ssf
    ssf_coefs <- apply(coef(fit_ssf)$cond$stratum, 2, mean)[-1]
    # Get variance-covariance matrix from ssf
    ssf_vcov <- vcov(fit_ssf)$cond
    # Sample from multivariate normal distribution of slopes/sigmas
    beta.hats <- MASS::mvrnorm(n = nsims, mu = ssf_coefs, Sigma = ssf_vcov)
    ntot.test <- nrow(xmat) # total number of test observations
    # Set up beta for generating prediction
    dimxmat <- ncol(xmat)
    tempdat <- as.data.frame(cbind(1:nrow(xmat),stratum, z))
    names(tempdat)[1] <- "inds"
    # Use new beta^ vector to generate predicted values for test data
    for(i in 1:nsims){
      lp.hat.s <- as.matrix(xmat)%*%beta.hats[i,]
      tempdat$wx.s <- exp(lp.hat.s)
      tempdat.g <- tempdat%>%group_by(stratum)
      x_sim_choice[i,,] <- z[sample_n(tempdat.g,1, weight=wx.s)$inds,]
    }
    return(x_sim_choice)
  }
  
  # Simulate from model x n_iterations
  xchoice <- simstrat(nsims = n_iterations,
                      xmat = design.mat.test, 
                      stratum = mdat.test$stratum, 
                      fit_ssf = ssf.train,
                      z = z)   
  # Get density estimates for habitat characteristics of observed locations in 
  # test data and associated with random locations generated by uhcsim
  # Initiate a df
  denshats_df <- data.frame()
  for(i in 1:(length(covs_no_ranef)+1)) {
    # Calculate density estimates
    denshats <- uhcdenscalc(rand_sims = xchoice[,,i], 
                            dat = z[mdat.test$presence==1,i], 
                            avail = z[mdat.test$presence==0,i],
                            gridsize = 500) 
    # Bind covariates with remaining data frame
    denscov <- data.frame(
      # Specify period
      period = calv_period,
      # Specify model type
      model_form,
      covariate = covs[i],
      # Get quantiles of the predicted distribution of the variable at used points
      densrand_l = apply(denshats[["densrand"]], 2,
                         function(x) quantile(x, probs = c(0.025, 0.975), na.rm = T))[1 ,],
      densrand_h = apply(denshats[["densrand"]], 2,
                         function(x) quantile(x, probs = c(0.025, 0.975), na.rm = T))[2 ,],
      # Get the actual variable distribution at used points
      densdat_x = denshats[["densdat"]]$x,
      densdat_y = denshats[["densdat"]]$y,
      # Get available distribution at used points
      densavail_x = denshats[["densavail"]]$x,
      densavail_y = denshats[["densavail"]]$y
    )
    # Bind together dataframes for all covariates
    denshats_df <- rbind(denshats_df, denscov)
  }
  
  return(denshats_df)
  
}


