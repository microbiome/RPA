\name{get.probe.noise.estimates}

\Rdversion{1.1}

\alias{get.probe.noise.estimates}

\title{Fetch probe-level noise estimates from an rpa object}

\description{Provides probe-level estimates of noise, as given by the
  RPA algorithm.}

\usage{get.probe.noise.estimates(rpa.res, sets = NULL, normalization = NULL, verbose = FALSE)}

\arguments{

  \item{rpa.res }{An rpa object.}

  \item{sets }{Probesets to check.}
  
  \item{normalization }{Normalization method for probe noise
    estimates.}

    The higher the value, the higher the probe-level noise. By default,
    the estimated probe-level variances (sigma2) of the RPA model are
    returned. Other options include:

    "withinset.weights": inverse of probe-wise weights used to summarize
    the probe-level observations into probeset-level signal
    
    "withinset.relative": probe-wise standard deviations versus standard
    deviation of the probeset-level signal d. This compensates for the
    effect that the detected probe-level noise is typically coupled with
    overall signal levels.

    "withinset.categorical": indexes the probes according to their
    reliability within each probeset. Probes with higher indices are
    more noisy.
    
  \item{verbose }{Print progress information during computation.}

}

\details{Provides probe-specific noise estimates. The normalization
  options are included to improve comparability across probesets.}

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

