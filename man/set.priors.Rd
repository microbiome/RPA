\name{set.priors}
\Rdversion{1.1}
\alias{set.priors}

\title{Set prior parameter object for RPA.}

\description{Reset some of the existing priors, or create a template of
priors for the whole data.}

\usage{set.priors(abatch, set, alpha, beta, priors = NULL, alpha.template = 1e-6, beta.template = 1e-6)}

\arguments{

\item{abatch }{An AffyBatch object.}
	    
\item{set }{Probeset where priors are to be determined.}

\item{alpha, beta }{Vectors giving the values for the probe-wise prior
parameters.}

\item{priors }{Optional. Existing prior to be modified.}

\item{alpha.template, beta.template }{Scalars. Can be used to define the
  default prior values for the whole-data prior template.}

}

\details{The method returns a prior object that specifies probe-wise
  priors for the whole data set. By default, it sets non-informative
  priors for all probes, except those specified by the parameters 'set',
  'alpha' and 'beta'. If a prior object is given in the input, then only
  the values for the specified probeset ('set') will be modified.}

\value{An instance of 'rpa.priors' class.}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE, to appear. See
http://www.cis.hut.fi/projects/mi/software/RPA/}

\author{Leo Lahti <leo.lahti@iki.fi>}

\examples{

require(affy)
require(affydata)
data(Dilution)

# Create a prior object with specific alpha, beta for one probeset
alpha <- beta <- rep(1, 16)
alpha[[5]] <- 3; beta[[5]] <- 1
priors <- set.priors(Dilution, set = "1000_at", alpha, beta)

## Run RPA using the predefined priors
# rpa.results <- RPA.pointestimate(Dilution, set, priors = priors)

}

\keyword{ methods }

