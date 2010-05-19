\name{get.probe.noise.estimates}

\Rdversion{1.1}

\alias{get.probe.noise.estimates}

\title{Fetch probe-level noise estimates from an rpa object}

\description{Provides estimates of probe-level noise, given by the
  RPA algorithm.}

\usage{get.probe.noise.estimates(rpa.res, sets = NULL, normalization = NULL, verbose = FALSE)}

\arguments{

  \item{rpa.res }{An rpa object.}

  \item{sets }{Probesets to check.}
  
  \item{normalization }{Optional normalization for probe noise
    estimates.}

    The higher the value, the higher the probe-level noise. By default,
    probe-level variances of the RPA model are returned. Other options
    include:

    "withinset.weights": The relative weight of a probe within probeset
       is determined by the relative noise of the probe with respect to
       the other probes in the same probeset. This option returns the
       inverse of probe-specific weights within each probeset. This can
       be used to normalize probe-level weights to improve comparability
       across probesets.

       "withinset.relative": The detected probe-level noise can be
         coupled with overall signal levels of the probeset. This option
         provides an estimate of probe-wise standard deviation
         normalized by the standard deviation of the probeset-level
         signal d.

	 "withinset.categorical": In some applications it can be
    	   sufficient to investigate the relative order of the probes,
    	   ignoring the parameter estimates. This option indexes the
    	   probes according to their reliability within each
    	   probeset. Probes with higher indices are more noisy.
    
  \item{verbose }{Print progress information during computation.}

}

\details{The normalization options are included to improve comparability
  across probesets. The higher the variance, the more noisy the
  probe. Inverse of the variance, can be used to quantitate probe
  reliability. Note that the relative weight of a probe within probeset
  is determined by the relative noise of the probe with respect to the
  other probes in the same probeset. Comparison of probe-specific
  variances across probesets may benefit from normalization of this
  effect. Therefore optional normalizations for probe noise estimation
  are provided.}

\value{A list. Each element corresponds to one probeset (of the input
  object). The element lists noise estimates for each probe within the
  probeset.}

\references{Probabilistic Analysis of Probe Reliability in
Differential Gene Expression Studies with Short Oligonucleotide
Arrays. Lahti et al., TCBB/IEEE. See
http://www.cis.hut.fi/projects/mi/software/RPA/ }

\author{Leo Lahti <leo.lahti@iki.fi>}

\seealso{RPA.pointestimate}

\examples{

## Load example data set
require(affydata)
data(Dilution)
## Compute RPA
rpa.results <- RPA.pointestimate(Dilution, set = "1000_at")
noise <- get.probe.noise.estimates(rpa.results)
}

\keyword{ methods }

