RPA.tau2.update <- function (R, alpha, beta, tau2.method = "robust") {

  # FIXME: online mode not necessarily needed at all here
  # FIXME: speedup by defining the method outside this looping function

  # tau2 values are notably smaller with 
  # tau2.method = "fast" than tau2.method = "robust"
  # This is because different priors for alpha are used
  # they lead to similar probe weights, though so not so much effect on probeset level estimates

  # R <- S - d: arrays x probes matrix; observations vs. estimated real signal
  # Note: alpha here is alphahat = T/2 + alpha w.r.t. user-defined alpha prior

  # alpha > 1 required in tau2.method 'mode' and 'robust'
  if (tau2.method == "mean" || tau2.method == "online" || tau2.method == "robust") {
    # mean for invgam(sig^2 | alpha,beta)
    tau2 <- beta / (alpha - 1) # FIXME: speedup by precalculating alpha - 1?
  } else if (tau2.method == "mode") {
    # mode for invgam(sig^2 | alpha,beta)
    tau2 <- beta / (alpha + 1) 
  } else if (tau2.method == "fast") {
    # this is often also used
    tau2 <- beta / alpha 
  } else if (tau2.method == "var") {
    # Assume uninformative priors alpha, beta -> 0	  
    # NOTE: RPA converges to variance with large sample size
    # priors not used
    # Do not center: S - d is already assumed to be centered by definition
    # S = d + N(0, sigmaj2)  
    # R <- S - d 
    # faster way to calculate variance than apply(R, 2, var):
    # Note that R is assumed to be centered already,
    # so actually this is more accurate for our case than apply(R,2,var)
    # which would center the data
    #tau2 <- colSums(centerData(R)^2)/(nrow(R) - 1) 
    tau2 <- colSums(R^2)/(nrow(R) - 1) 
  } else { 
    stop("Invalid tau2.method provided!") 
  }

  # return updated tau2
  tau2
}


