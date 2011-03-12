\name{RPA-package}

\Rdversion{1.1}

\alias{RPA-package}

\alias{RPA}

\docType{package}

\title{RPA: probabilistic analysis of probe reliability and gene expression}

\description{RPA estimates probe-specific variances and differential
gene expression using probe-level observations of differential gene
expression.}

\details{  
  \tabular{ll}{
Package: \tab RPA\cr
Type: \tab Package\cr
Version: \tab 1.7.32\cr
Date: \tab 2011-03-12\cr
License: \tab FreeBSD\cr
LazyLoad: \tab yes\cr
}
RPA.pointestimate is used to calculate probe reliability and
differential expression estimates: 'rpa.results <-
RPA.pointestimate(affybatch)'. The other functions are provided for
users who wish to investigate the details of the algorithm more
closely.
}
\author{Leo Lahti \email{leo.lahti@iki.fi}}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE 2011.}

\keyword{ package }

\examples{

  ## Not run:

  ## Load example data set (Dilution affybatch).
  ## This is a toy example with a small example dataset
  ## for probe reliability analysis (4 arrays).
  ## For practical applications, a larger sample size is
  ## recommended.

  #require(affy)
  #require(affydata)
  #data(Dilution)

  ## Compute RPA for whole data set
  ## ... slow, not executed here
  #rpa.results <- RPA.pointestimate(Dilution)

}
