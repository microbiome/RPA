
set.alpha <- function (alpha, sigma2.method, P){ 

  # if alpha is scalar, set identical prior for all probes with this value
  if (is.null(alpha) && !(sigma2.method == "robust" || sigma2.method == "mode" || sigma2.method == "online")) {
    # uninformative
    alpha <- rep.int(1e-6, P)
  } else if (is.null(alpha) && (sigma2.method == "robust" || sigma2.method == "robust")) {
    # alpha not given: set equal and informative priors to
    # avoid collapse to individual probes
    alpha <- rep.int(2, P)
  } else if (is.null(alpha) && sigma2.method == "mean") {
    # alpha not given: set equal and informative priors to
    # avoid collapse to individual probes
    alpha <- rep.int(1 + 1e-6, P) # alpha > 1 required
  } else if (length(alpha) == 1) {
    alpha <- rep.int(alpha, P)
  } else {
    alpha <- alpha
  }

  if ((sigma2.method == "mean" || sigma2.method == "online" || sigma2.method == "robust") && any(alpha <= 1)) {
    alpha <- rep(1+1e-6, P)
    warning(paste("Initial alpha parameter has to be >1. Setting alpha to:", alpha))
  }

  alpha
}

set.beta <- function (beta, sigma2.method, P) {

  # if beta is scalar, set identical prior for all probes with this value
  if (is.null(beta) && !sigma2.method == "robust") {
    beta <- rep.int(1e-6, P)
  } else if (is.null(beta) && sigma2.method == "robust") {
    beta <- rep.int(1, P)
  } else if (length(beta) == 1) {
    beta <- rep.int(beta, P)
  } else {}

  beta

}

##################################

centerData <- function (X,rm.na = FALSE, meanvalue = NULL) {


  if (!rm.na) {
    xcenter <- colMeans(X)
    X2 <- X - rep(xcenter, rep.int(nrow(X), ncol(X)))
  } else {	
    X2 <- array(NA, dim = c(nrow(X), ncol(X)), dimnames = dimnames(X))
    for (i in 1:ncol(X)) {
      x <- X[,i]
      nainds <- is.na(x)
      xmean <- mean(x[!nainds])
      X2[!nainds,i] <- x[!nainds] - xmean 	
    }
    dimnames(X2) <- dimnames(X)
  }

  if (!is.null(meanvalue)) {
    # Shift the data so that mean gets a specified value
    X2 <- X2 + meanvalue
  }


  
  X2
}

