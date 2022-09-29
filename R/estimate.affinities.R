#' @title estimate.affinities
#' @description Probe affinity estimation. Estimates probe-specific affinity parameters.
#' 
#' @param dat Input data set: probes x samples.
#' @param a Estimated expression signal from RPA model.
#'
#' @export
#'
#' @details 
#'  To estimate means in the original data domain let us assume that
#'  each probe-level observation x is of the following form:
#'  x = d + v + noise,
#'  where x and d are vectors over samples,
#'  v is a scalar (vector with identical elements)
#'  noise is Gaussian with zero mean and probe-specific variance parameters tau2 
#'  Then the parameter mu will indicate how much probe-level observation
#'  deviates from the estimated signal shape d.
#'  This deviation is further decomposed as
#'  mu = mu.real + mu.probe, where
#'  mu.real describes the 'real' signal level, common for all probes
#'  mu.probe describes probe affinity effect
#'  Let us now assume that mu.probe ~ N(0, sigma.probe).
#'  This encodes the assumption that in general the affinity
#'  effect of each probe tends to be close to zero.
#'  Then we just calculate ML estimates of mu.real and mu.probe
#'  based on particular assumptions.
#'  Note that this part of the algorithm has not been defined
#'  in full probabilistic terms yet, just calculating the point estimates.

#'  Note that while tau2 in RPA measures stochastic noise, and NOT the
#'  affinity effect, we use it here as a heuristic solution to weigh the
#'  probes according to how much they contribute to the overall signal
#'  shape. Intuitively, probes that have little effect on the signal
#'  shape (i.e. are very noisy and likely to be contaminated by many
#'  unrelated signals) should also contribute less to the absolute
#'  signal estimate. If no other prior information is available, using
#'  stochastic parameters tau2 to determine probe weights is likely to
#'  work better than simple averaging of the probes without
#'  weights. Also in this case the probe affinities sum close to zero
#'  but there is some flexibility, and more noisy probes can be
#'  downweighted.
#'
#' @return A vector with probe-specific affinities.
#'
#' @seealso rpa.fit
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples #  mu <- estimate.affinities(dat, a)
#' @keywords utilities

estimate.affinities <- function (dat, a) {

  # In RPA, the final signal estimate
  # is weighted average over the probes, 
  # both in terms of shape and affinities
  
  # Calculate affinities w.r.t. overall signal
  colMeans(t(dat) - a)

}
