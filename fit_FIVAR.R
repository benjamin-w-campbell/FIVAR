######################################
# fit_FIVAR.R
# Aut: BW Campbell
# Date: 01 Dec. 2016
# Purpose: Estimate a
#   Fractionally Integrated
#   VAR estimated according
#   to the regime outlined
#   in Box-Steffensmeier et. al (2009)
######################################

# This function follows a fairly clear set of steps
# First: It loops through the endogenous variables and pre-whitens them according to
  # the best fitting ARFIMA model, storing the residuals once all ARMA(p,q) are accounted for.
  # It then differences those residuals according to the degree of fractional integration estimated.
# Second, it selects the optimal number of lags according to model AIC.
# Third, a VAR model is estimated with the optimal number of lags on the pre-whitened data.

# Feel free to modify the function as you see fit; you may decide to select your own lags instead 
  # of having them automatically selcted according to AIC.    

# fit_FIVAR takes a few different inputs
  # endog_vars: the variables that are to be specified as endogenous in the VAR model. 
  # drange: the range that you allow the fractional integration parameter "d" to take on.  Limiting ensures model stability.
  # max_lags: as lags are to be selected by model fit, one must determine the maximum number of lags that the VAR selection
    # process can take.  Depending upon the number of time steps, this may be more or less than 10.
  # exog_vars: exogeneous variables for the VAR model.

fit_FIVAR <- function(endog_vars, drange = c(0, 0.85), max_lags = 10, exog_vars = NULL){
  # load relevant packages
  require(forecast)
  require(fracdiff)
  require(vars)
  
  # put the ts_dat object, which may be a concatenated vector as a list
  dat <- as.list(endog_vars)
  
  # Pre-whiten the data:
  # if this list is of length one, it changes slightly the return process
  if(length(dat) == 1){
    d_out <- list()
    endog_vars <- list()
    fit <- arfima(dat[[1]], drange = drange)
    pre_whitening <- diffseries(fit$residuals, d = fit$d)
    endog_vars <- pre_whitening
    d_out  <- fit$d
  }
  
  if(length(dat) > 1){
    d_out <- list()
    endog_vars <- list()
    for(i in 1:length(dat)){
      fit <- arfima(dat[[i]], drange = drange)
      pre_whitening <- diffseries(fit$residuals, d = fit$d)
      endog_vars[[i]] <- pre_whitening
      d_out[[i]] <- fit$d
    }
  }
  
  endog_vars_df <- as.data.frame(unlist(endog_vars[[1]]))
  colnames(endog_vars_df)[1] <- "X1"
  for(i in 2:length(endog_vars)){
    vec <- as.data.frame(unlist(endog_vars[[i]]))
    endog_vars_df <- cbind(endog_vars_df, vec)
    colname <- paste0("X", i)
    colnames(endog_vars_df)[i] <- colname
  }
  
  # Now, fit the VAR to the pre-whitenned series
  selectedLags <- VARselect(y = endog_vars_df,
                            lag.max = max_lags,
                            exogen = exog_vars)
  
  # AIC is used to select lags
  lags <- unlist(selectedLags$selection[1])
  
  # Estimate VAR
  mod <- VAR(endog_vars_df, p = lags, type = "const", exog = exog_vars)
  
  detach("package:forecast", unload = TRUE)
  detach("package:fracdiff", unload = TRUE)
  detach("package:vars", unload = TRUE)
  out_list <- list(endog_vars, lags, d_out, selectedLags, mod)
  
 return(out_list)
}

