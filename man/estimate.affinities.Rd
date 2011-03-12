\name{estimate.affinities}
\alias{estimate.affinities}
\title{Probe affinity estimation}
\description{Estimates probe-specific affinity parameters.}
\usage{
estimate.affinities(dat, d, sigmas = NULL)
}

\arguments{
  \item{dat}{Input data set: probes x samples.}
  \item{d}{Estimated differential expression signal from RPA
  model. Determines the shape of the true signal in the model.}
  \item{sigmas}{Probe-specific priors for affinity parameters in terms
  of N(0, sigma2.affinity) model. If not given, identical prior is given
  for all probes.}
}
\details{

  To estimate means in the original data domain let us assume that
  each probe-level observation x is of the following form:
  x = d + mu + noise,
  where x and d are vectors over samples,
  mu is a scalar (vector with identical elements)
  noise is Gaussian with zero mean and probe-specific variance parameters sigma2 
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

  if identical sigma.probe is used for all probes then mu.real is
  estimated by the average of the probe effects mu. Note that then
  probe-specific affinities mu.probe will sum to exactly zero, which
  gives an analogous model than used in RMA, which uses this
  assumption to fit medianpolish.

  Probes can also be weighted by setting probe-specific sigmas.  In
  "rpa" option, set sigma.probe to the sigma2 value of the probe
  estimated by RPA.  Note that while sigma2 in RPA measures
  stochastic noise, and NOT the affinity effect, we can use it here
  as a heuristic solution to weigh the probes according to how much
  they contribute to the overall signal shape. Intuitively, probes
  that have little effect on the signal shape (i.e. are very noisy
  and likely to be contaminated by many unrelated signals) should
  also contribute less to the absolute signal estimate. If no other
  prior information is available, using stochastic parameters sigma2
  to determine probe weights is likely to work better than simple
  averaging of the probes without weights. Also in this case the
  probe affinities sum close to zero but there is some flexibility,
  and more noisy probes can be downweighted.

}
\value{
  A list with the following elements:
  \item{real }{Estimated 'real' signal level (on top of differential
    expression signal).}
  \item{probe }{Probe-specific affinities. A list with each element
    corresponding to one probe.}
}


\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE. See
http://bioconductor.org/packages/release/bioc/html/RPA.html}

\author{Leo Lahti \email{leo.lahti@iki.fi}}



\note{Affinity estimation is not part of the original RPA procedure in
TCBB/IEEE 2011 paper. It is added here since estimates of the absolute
levels are often needed in microarray applications. Note that affinity
parameters are unidentifiable in the model if no prior assumptions are
given. We assume that affinity effects are zero on average, but allow
some flexibility through probe-specific weights.}

\seealso{rpa.fit}
\examples{
##  mu <- estimate.affinities(dat, d, sigmas = NULL)
}
\keyword{ utilities }

