\name{RPA.dcost}
\Rdversion{1.1}
\alias{RPA.dcost}
\title{Cost function for d.}
\description{Internal function.}
\usage{RPA.dcost(d, sigma2, S)}

\arguments{
  
  \item{d}{A vector. Estimated 'true' signal underlying the noisy
  probe-level observations.}
  
  \item{sigma2}{A vector of probe-specific variances.}

  \item{S}{Matrix of probe-level observations for a single probeset:
    samples x probes.}
}

\value{Cost for d given S and sigma2.}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE, to appear. See
http://www.cis.hut.fi/projects/mi/software/RPA/}

\author{Leo Lahti <leo.lahti@tkk.fi>}

\seealso{RPA.iteration}
\keyword{ internal }

