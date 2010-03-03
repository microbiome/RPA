\name{RPA-package}

\Rdversion{1.1}

\alias{RPA-package}

\alias{RPA}

\docType{package}

\title{RPA: probe reliability and differential expression analysis}

\description{RPA estimates probe-specific variances and differential
gene expression using probe-level observations of differential gene
expression.}

\details{  
  \tabular{ll}{
Package: \tab RPA\cr
Type: \tab Package\cr
Version: \tab 1.3.2\cr
Date: \tab 2010-03-02\cr
License: \tab GNU GPL 2 or any later version (at your option)\cr
LazyLoad: \tab yes\cr
}
RPA.pointestimate is used to calculate probe reliability and
differential expression estimates: 'rpa.results <-
RPA.pointestimate(affybatch)'. The other functions are provided for
users who wish to investigate the details of the algorithm more
closely.
}
\author{Leo Lahti <leo.lahti@iki.fi>}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE. See
http://www.cis.hut.fi/projects/mi/software/RPA/ }


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
