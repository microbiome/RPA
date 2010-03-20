get.probe.noise.estimates <- function (rpa.res, sets = NULL, normalization = NULL, verbose = FALSE) {

  if (is.null(sets)) {sets <- names(rpa.res)}
  
  if ( is.null( normalization ) ) {
    # No normalization: return probe variances
    rpa.res$sigma2 <- rpa.res$sigma2

  } else if (normalization == "withinset.categorical") {

    for (set in sets) {
      # Order probes within the probeset by their noise level
      # smallest index means least noise
      if (verbose) {cat(set)}
      o <- order(rpa.res$sigma2[[set]])
      rpa.res$sigma2[[set]][o] <- 1:length(o)      
    }
    
  } else if (normalization == "withinset.relative") {

    for (set in sets) {
      # compare probe standard deviation to probeset signal standard
      # deviation
      # within the probeset
      # smallest means least noise
      if (verbose) {cat(set)}
      vars <- rpa.res$sigma2[[set]]
      rpa.res$sigma2[[set]] <- sqrt(vars)/sqrt(var(rpa.res$d[set,]))
    }
  } else if (normalization == "withinset.weights") {

    for (set in sets) {
      # inverse of the weights used to compute
      # probeset-level summary signal d
      # within the probeset
      # smallest means least noise (largest weight = smallest inverse weight)
      if (verbose) {cat(set)}
      vars <- rpa.res$sigma2[[set]]
      rpa.res$sigma2[[set]] <- 1/((1 / vars) / sum(1 / vars))
    }
  }

  rpa.res$sigma2[sets]

}
