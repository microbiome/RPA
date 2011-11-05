get.set.inds <- function (cel.files, cdf, sets) {

  # Get probe position indices
  abatch <- ReadAffy(filenames = cel.files[1:2], compress=getOption("BioC")$affy$compress.cel) 
  # Set alternative CDF environment if given
  if (!is.null(cdf)) { abatch@cdfName <- cdf }    

  # Check names for the investigated probesets
  # if my.sets not specified, take all sets in abatch
  if ( is.null(sets) ) { sets <- geneNames(abatch) } 

  # Retrieve probe positions
  # The indices correspond to the output of pm()
  pN <- probeNames(abatch)
  set.inds <- split(1:length(pN), pN)[sets] # pNList
 
  list(set.inds = set.inds, cdf = cdf)
 
}



get.batches <- function (items, batch.size, shuffle = FALSE) {

  # Random ordering for the items?
  if (shuffle) {items <- sample(items)}

  # N elements into batches of size batch.size
  # last batch can be smaller

  if (length(items) == 1 && is.numeric(items)) {
    N <- items
    items <- 1:N
  } else {
    N <- length(items)
  }

  if (N < batch.size) {
     warning("batch.size > N, setting batch.size = N.")
  }

  batches <- list()
  ns <- floor(c(0, seq(batch.size, N, batch.size)))
  if (ns[[length(ns)]] < N) {ns[[length(ns) + 1]] <- N}
  cnt <- 0
  for (i in 2:length(ns)) {
    cnt <- cnt + 1
    n.start <- ns[[i-1]]
    n.stop <- ns[[i]]
    batches[[i-1]] <- items[(n.start+1):n.stop]
  }
  
  batches
}


set.alpha <- function (alpha, sigma2.method, P){ 

  # if alpha is scalar, set identical prior for all probes with this value
  if (is.null(alpha)) {
    if (sigma2.method == "mode" || sigma2.method == "var") {
      alpha <- 1e-6     # uninformative
    } else if (sigma2.method == "robust" || sigma2.method == "mean" || sigma2.method == "online") {
      # alpha not given: set equal and informative priors to
      # avoid collapse to individual probes
      alpha <- 2
    } 
  }

  if ((sigma2.method == "mean" || sigma2.method == "online" || sigma2.method == "robust") && any(alpha <= 1)) {
    stop(paste("Set alpha > 1!"))
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

