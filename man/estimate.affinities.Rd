\name{estimate.affinities}
\alias{estimate.affinities}
\title{Probe affinity estimation}
\description{Estimates probe-specific affinity parameters.}
\usage{
estimate.affinities(dat, a)
}

\arguments{
  \item{dat}{Input data set: probes x samples.}
  \item{a}{Estimated expression signal from RPA model.}
}
\details{

  To estimate means in the original data domain let us assume that
  each probe-level observation x is of the following form:
  x = d + v + noise,
  where x and d are vectors over samples,
  v is a scalar (vector with identical elements)
  noise is Gaussian with zero mean and probe-specific variance parameters tau2 
  Then the parameter mu will indicate how much probe-level observation
  deviates from the estimated signal shape d.
  This deviation is further decomposed as
  mu = mu.real + mu.probe, where
  mu.real describes the 'real' signal level, common for all probes
  mu.probe describes probe affinity effect
  Let us now assume that mu.probe ~ N(0, sigma.probe).
  This encodes the assumption that in general the affinity
  effect of each probe tends to be close to zero.
  Then we just calculate ML estimates of mu.real and mu.probe
  based on particular assumptions.
  Note that this part of the algorithm has not been defined
  in full probabilistic terms yet, just calculating the point estimates.

  Note that while tau2 in RPA measures stochastic noise, and NOT the
  affinity effect, we use it here as a heuristic solution to weigh the
  probes according to how much they contribute to the overall signal
  shape. Intuitively, probes that have little effect on the signal
  shape (i.e. are very noisy and likely to be contaminated by many
  unrelated signals) should also contribute less to the absolute
  signal estimate. If no other prior information is available, using
  stochastic parameters tau2 to determine probe weights is likely to
  work better than simple averaging of the probes without
  weights. Also in this case the probe affinities sum close to zero
  but there is some flexibility, and more noisy probes can be
  downweighted.

}
\value{A vector with probe-specific affinities.}


\references{See citation("RPA").}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\note{Affinity estimation is not part of the original RPA procedure in
TCBB/IEEE 2011 paper. It is added here since estimates of the absolute
levels are often needed in microarray applications. Note that affinity
parameters are unidentifiable in the model if no prior assumptions are
given. We assume that affinity effects are zero on average, but allow
flexibility through probe-specific weights.}

\seealso{rpa.fit}
\examples{
##  mu <- estimate.affinities(dat, a)
}
\keyword{ utilities }

