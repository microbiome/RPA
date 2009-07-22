\name{d.update.fast}
\Rdversion{1.1}
\alias{d.update.fast}
\title{Fast update for probeset-level signal d.}

\description{Weighted average over the columns of S, weighted by
the inverse variances 1/sigma2.}

\usage{d.update.fast(S, sigma2)}

\arguments{
  
  \item{S}{Matrix of probe-level observations for a single probeset: samples x probes.}

  \item{sigma2}{A vector. Estimated variance for each measurement (probe).}
  
}

\value{A vector giving the estimated true signal underlying the
observations in S. }

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE, to appear. See
http://www.cis.hut.fi/projects/mi/software/RPA/}

\author{Leo Lahti <leo.lahti@tkk.fi>}

\examples{
set.seed(24)
vec <- seq(10)
S <- cbind(rnorm(length(vec),sd=1),
	  rnorm(length(vec),sd=2),
	  rnorm(length(vec),sd=3))
d.update.fast(S,sigma2=seq(3))
}
\keyword{ methods }
\keyword{ iteration }

