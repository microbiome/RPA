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
Version: \tab 1.1.2\cr
Date: \tab 2009-07-22\cr
License: \tab GNU GPL 2 or any later version (at your option)\cr
LazyLoad: \tab yes\cr
}
RPA.pointestimate computes probe reliability and differential expression
estimates: 'rpa.results <- RPA.pointestimate(affybatch)'. The other
functions are provided for users who wish to investigate the details of
the algorithm more closely.
}
\author{Leo Lahti <leo.lahti@tkk.fi>}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE, to appear. See
http://www.cis.hut.fi/projects/mi/software/RPA/ }


\keyword{ package }

\examples{
        ## Load example data set (Dilution affybatch).
        ## This is a toy example with a small example dataset
        ## for probe reliability analysis (4 arrays).
	## For practical applications, a larger sample size is
        ## recommended.

        require(affy)
        require(affydata)
        data(Dilution)

        ## Run RPA analysis
	## Compute RPA for the whole data set
	## Slow, not executed here
	##rpa.results <- RPA.pointestimate(Dilution)

	# Compute RPA for specific probesets only
	sets = geneNames(Dilution)[1:5]
	rpa.results <- RPA.pointestimate(Dilution,sets)

	## Visualize the results for one of the probe sets
	plot(rpa.results[sets[[1]],])
	
}
