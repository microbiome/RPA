#' @title rpa.summarize
#' @description RPA summarization.
#' @param dat Original data: probes x samples.
#' @param affinities Probe affinities
#' @param variances Probe variances
#' @param summarize.with.affinities Use affinity estimates in probe summarization step. Default: FALSE.
#' @details Summarizes the probes in a probe set according to the RPA model based on the given affinity and variance parameters.
#'
#' @return A vector. Probeset-level summary signal.
#'
#' @seealso rpa
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # res <- rpa.summarize(dat, affinities, variances, summarize.with.affinities = FALSE)
#' @keywords utilities

rpa.summarize <- function (dat, affinities, variances, summarize.with.affinities = FALSE) {

  # FIXME: add hyperparameter estimation for the case where 
  # affinities & variances are NULL

  # Impute if there are missing values
  dat <- rpa.impute(dat)
  if ( is.null(colnames(dat)) ) { colnames(dat) <- 1:ncol(dat) }

  # Accommodate single-probe probesets
  if (nrow(dat) == 1) {  

    mu <- as.vector(dat)
    mu <- mu - affinities
    names(mu) <- colnames(dat)

  } else {

    # Remove affinities from raw signal before summarization

    # Since probe-specific variance is now known (from the estimation above), 
    # the probeset-level signal
    # estimate is obtained as a weighted sum of the 
    # probes, weighted by the probe-specific variances
    
    if (summarize.with.affinities) {
      mu <- d.update.fast(dat - affinities, variances)
    } else {
      mu <- d.update.fast(dat, variances) # ignore affinities in summarization    
    }
  }
  
  mu

}
