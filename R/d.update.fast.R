#' @title Fast d update
#' @description Computes weighted average over the probes, weighted by their inverse probe-specific variances.
#' @param St probes x samples data matrix
#' @param s2 variances for the probes
#'
#' @details Returns summarized probeset-level weighted average
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords methods

d.update.fast <- function (St, s2) {

  # With large sample sizes when T -> Inf
  # d converges to the weighted mean 
  # over the probes, weighted by probe variances	 
  if (nrow(St) == 1) {
    v <- St[1,]
  } else {
    v <- colSums(St / s2) / sum(1 / s2)
  }

  v

}


