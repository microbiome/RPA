#' @title RPA iteration
#' @description Estimating model parameters d and tau2.
#' @param S Matrix of probe-level observations for a single probeset: samples x probes.
#' @param epsilon Convergence tolerance. The iteration is deemed converged when the change in all parameters is < epsilon.
#' @param alpha alpha prior for inverse Gamma distribution of probe-specific variances. Noninformative prior is obtained with alpha, beta -> 0.  Not used with tau2.method 'var'. Scalar alpha and beta are specify equal inverse Gamma prior for all probes to regularize the solution. The defaults depend on the method.
#' @param beta beta prior for inverse Gamma distribution of probe-specific variances. Noninformative prior is obtained with alpha, beta -> 0.  Not used with tau2.method 'var'. Scalar alpha and beta are specify equal inverse Gamma prior for all probes to regularize the solution. The defaults depend on the method.
#' @param tau2.method Optimization method for tau2 (probe-specific variances).
#'
#'     "robust": (default) update tau2 by posterior mean,
#'		regularized by informative priors that are identical
#'		for all probes (user-specified by
#'		setting scalar values for alpha, beta). This
#'		regularizes the solution, and avoids overfitting where
#'		a single probe obtains infinite reliability. This is a
#'	        potential problem in the other tau2 update
#'	        methods with non-informative variance priors. The
#'		default values alpha = 2; beta = 1 are
#'	        used if alpha and beta are not specified.
#' 
#'      "mode": update tau2 with posterior mean
#'
#'	"mean": update tau2 with posterior mean
#'	
#'	"var": update tau2 with variance around d. Applies the fact
#'             that tau2 cost function converges to variance with
#'               large sample sizes. 
#'
#' @param d.method Method to optimize d.
#'        "fast": (default) weighted mean over the probes, weighted by
#'		probe variances The solution converges to this with
#'		large sample size.
#'
#'        "basic": optimization scheme to find a mode used in Lahti et
#'        	 al. TCBB/IEEE; relatively slow; this is the preferred 
#'		 method with small sample sizes.
#'                	      
#' @param maxloop Maximum number of iterations in the estimation process.
#'
#' @details Finds point estimates of the model parameters d (estimated true signal underlying probe-level observations), and tau2 (probe-specific variances). Assuming data set S with P observations of signal d with Gaussian noise that is specific for each observation (specified by a vector tau2 of length P), this method gives a point estimate of d and tau2. Probe-level variance priors alpha, beta can be used with tau2.methods 'robust', 'mode', and 'mean'.  The d.method = "fast" is the recommended method for point computing point estimates with large samples size.
#'
#' @return A list with the following elements: d:  A vector. Estimated 'true' signal underlying the noisy probe-level observations.; tau2: A vector. Estimated variances for each measurement (or probe).
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords utilities

RPA.iteration <- function(S,
                          epsilon = 1e-3,
                            alpha = NULL,
                             beta = NULL,
                    tau2.method = "fast",
                         d.method = "fast",
                          maxloop = 1e6)
{


  # FIXME: remove d.method and tau2.method from options
  # S <- q; maxloop <- 1e6

  P <- ncol(S) # number of probes
  T <- nrow(S) # Number of arrays (except reference)

  # Accommodate single-probe probesets
  if (P == 1) { return(list(d = as.vector(S), tau2 = 0, alpha = T/2, beta = 0)) }

  # Check: if affybatch/probeset is erroneous and contains just NAs or NaNs then return NA vector
  if (all(is.nan(S) | is.na(S))) { 
    return(list(d = rep(NA, T), tau2 = rep(NA, P))) 
  }

  # uninformative priors for tau2.methods mean, mode, var;
  # informative for 'robust', or alpha, beta are provided by user
  alpha.prior <- alpha <- set.alpha(alpha, tau2.method, P)
  beta.prior <- beta <- set.beta(beta, tau2.method, P)

  # Confirm that alpha is valid for tau2.method 
  if (tau2.method == "mean" || tau2.method == "robust") {
    ifelse(all(alpha > 1), TRUE, stop("alpha > 1 - (N.arrays - 1) / 2 required for this tau2.method"))
  } 

  ###############################

  # initialize tau2 with user-defined priors

  if (tau2.method == "var") {
    s2.meth <- "mean"
  } else {
    s2.meth <- tau2.method
  }

  # FIXME: add chance to give user-defined priors in here!
  tau2 <- RPA.tau2.update(NULL, alpha.prior, beta.prior, s2.meth)
  # initialize with large initial variance if not 
  if (length(tau2) == 0) {tau2 <- rep(3*sd(as.vector(S)), P)} 

  # check convergence at first iteration, 
  # not used in calculations
  tau2.old <- rep(Inf, P)
  
  ###############################

  # Update alpha
  # Do NOT update beta yet; it will be updated in while loop
  # after estimating d based on the current priors!
  alpha <- update.alpha(T, alpha) 

  #################################

  # optimize until convergence
  loopcnt <- 0

  # initialize d for the first iteration
  d <- rowMeans(S)  

  while ((max(abs(c(tau2 - tau2.old))) > epsilon) && loopcnt < maxloop) {

      tau2.old <- tau2

      # update d, given tau2
      d <- d.update.fast(t(S), tau2) # d.method = "fast"
      # d <- optim(d, fn = RPA.dcost, method = "BFGS", tau2 = tau2, S = S)$par # d.method = "basic"

      # Estimate noise 
      R <- S - d

      # beta update (feed in beta prior, no updates from this loop!)
      beta <- update.beta(R, beta.prior)

      # update tau2
      tau2 <- RPA.tau2.update(R, alpha, beta, tau2.method)

      # follow iteration count to avoid potentially infinite loops
      loopcnt <- loopcnt + 1 
  }

  list(d = d, tau2 = tau2, alpha = alpha, beta = beta)
}

#################################################################################

# Changelog:

RPA.iteration.fast <- function(S,
                          epsilon = 1e-3,
                            alpha.prior,
			    alpha.posterior,
                             beta.prior,
                          maxloop = 1e6)
{

  # S: samples x probes 

  # Check: if affybatch/probeset is erroneous and contains just NAs or NaNs then return NA vector
  if (all(is.nan(S) | is.na(S))) { 
    return(list(d = rep(NA, nrow(S)), tau2 = rep(NA, ncol(S)))) 
  }

  # Initialize variances by the original user-specified priors 
  tau2 <- beta.prior/alpha.prior # approx. mean and mode when sample size is largish

  # check convergence at first iteration, 
  # not used in calculations
  tau2.old <- rep(-1e6, length(tau2))
  
  ###############################

  # optimize until convergence
  loopcnt <- 0

  while ((max(abs(c(tau2 - tau2.old))) > epsilon) && loopcnt < maxloop) {

      tau2.old <- tau2

      # update d, given tau2
      d <- d.update.fast(t(S), tau2)

      # Estimate noise 
      R <- S - d

      # beta update (feed in original beta prior, not updates from this loop!)
      beta.posterior <- beta.fast(beta.prior, R) 
      #betahat.f(beta.prior, R) # update.beta(R, beta.prior)

      # update tau2 / NOTE: use posterior alpha and beta here!
      tau2 <- beta.posterior/alpha.posterior
      #RPA.tau2.update(R, alpha.posterior, beta.posterior, tau2.method = fast)

      # follow iteration count to avoid potentially infinite loops
      loopcnt <- loopcnt + 1 

  }


  list(d = d, tau2 = tau2, beta = beta.posterior)
}
