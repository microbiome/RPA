library(RPA)

fs <- list.files("../RPA/R", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}

require(affy)
require(affydata)
data(Dilution)

#eset <- rpa(Dilution)
#require(RPA)
sets <- geneNames(Dilution)[1:2]

#rpa.results <- RPA.pointestimate(Dilution, sets)
abatch<-Dilution; myseed = 101; priors = NULL; epsilon = 1e-2; cind = 1; sigma2.method = "robust"; d.method = "fast"; verbose = TRUE; bg.method = "rma"; normalization.method = "quantiles.robust"; cdf = NULL

  # PREPROCESSING

  #################################################################

  #Set random seed
  set.seed( myseed )

  # Set alternative CDF environment if given
  if (!is.null(cdf)) { abatch@cdfName <- cdf }
      
  # Preprocessing
  preproc <- RPA.preprocess(abatch, bg.method, normalization.method, cdf)

  # Number of arrays 
  T <- ncol(exprs(abatch))

  # Check names and number for the investigated probesets
  # if my.sets not specified, take all sets in abatch
  if ( is.null(sets) ) { sets <- geneNames(abatch) } 

  if (!all(sets %in% probeNames(abatch))) {
    warning("some probesets (sets) not available in abatch!")
    sets <- sets[sets %in% probeNames(abatch)]
  } else {}
  
  Nsets <- length(sets) 

  ## Matrices to store the results (also including reference array)
  d.results <- array(NA, dim = c(Nsets, T))
  rownames(d.results) <- sets
  colnames(d.results) <- colnames(exprs(abatch))

  sigma2.results <- vector(length = Nsets, mode = "list")  
  names(sigma2.results) <- sets

  affinity.results <- vector(length = Nsets, mode = "list")  
  names(affinity.results) <- sets

  mu.real <- vector(length = Nsets, mode = "list")  
  names(mu.real) <- sets    

  # Pick the priors for this set (gives NULL if no prior has been defined)
  alpha <- priors$alpha

  for (i in 1:Nsets) {

    set <- sets[[i]]
  
    if (verbose) {message(paste("Summarizing probeset", set, ":", i, "/", Nsets, "...\n"))}

    # Find probe (pm) indices for this set
    pmindices <- preproc$set.inds[[set]]

    # Pick the priors for this set (gives NULL if no prior has been defined)
    beta  <- priors[[set]]$beta
        
dat <- matrix(preproc$q[pmindices,], length(pmindices)); epsilon <- 1e-2
#    res <- rpa.fit(
#rpa.fit <- function (dat, cind = 1, epsilon = 1e-2, alpha = NULL, beta = NULL, sigma2.method = "fast", d.method = "fast") {

  if (is.null(colnames(dat))) {colnames(dat) <- 1:ncol(dat)}

  # Fit RPA
  S <- t(dat[, -cind] - dat[, cind])
#  estimated <- RPA.iteration(S, epsilon, alpha, beta, sigma2.method, d.method)

  P <- ncol(S) # number of probes
  T <- nrow(S) # Number of arrays (except reference)

  # Check: if affybatch/probeset is erroneous and contains just NAs or NaNs then return NA vector
  if (all(is.nan(S) | is.na(S))) { 
    return(list(d = rep(NA, T), sigma2 = rep(NA, P))) 
  }

  # uninformative priors for sigma2.methods mean, mode, var;
  # informative for 'robust', or alpha, beta are provided by user
  alpha.prior <- alpha <- set.alpha(alpha, sigma2.method, P)
  beta.prior <- beta <- set.beta(beta, sigma2.method, P)

  # Confirm that alpha is valid for sigma2.method 
  if (sigma2.method == "mean" || sigma2.method == "robust") {
    ifelse(all(alpha > 1), TRUE, stop("alpha > 1 - (N.arrays - 1) / 2 required for this sigma2.method"))
  } else {}

  ###############################

  # initialize sigma2 with user-defined priors

  if (sigma2.method == "var") {
    s2.meth <- "mean"
  } else {
    s2.meth <- sigma2.method
  }
  sigma2 <- RPA.sigma2.update(NULL, alpha.prior, beta.prior, s2.meth)

  # check convergence at first iteration, 
  # not used in calculations
  sigma2.old <- rep(-Inf, length(sigma2))
  
  ###############################

  # Update alpha
  # Do NOT update beta yet; it will be updated in while loop
  # after estimating d based on the current priors!
  alpha <- update.alpha(T, alpha) 

  #################################

  # optimize until convergence
  loopcnt <- 0

    while ((max(abs(c(sigma2 - sigma2.old))) > epsilon) && loopcnt < maxloop) {

      sigma2.old <- sigma2

      # update d, given sigma2
      d <- d.update.fast(t(S), sigma2)

      # Estimate noise 
      R <- S - d

      # beta update (feed in beta prior, not updates from this loop!)
      beta <- update.beta(R, beta.prior)

      # update sigma2
      sigma2 <- RPA.sigma2.update(R, alpha, beta, sigma2.method)

      # follow iteration count to avoid potentially infinite loops
      loopcnt <- loopcnt + 1 

    }

  
