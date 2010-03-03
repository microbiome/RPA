\name{rpa2eset}
\Rdversion{1.1}
\alias{rpa2eset}

\title{Coerce 'rpa' object into an 'ExpressionSet'}

\description{An instance of 'rpa' class contains differential gene
expression estimates in the variable 'd'. The function 'rpa2eset'
coerces this into an ExpressionSet object to allow downstream analysis
of the results using standard R/BioC tools for gene expression data.}

\usage{rpa2eset(x)}

\arguments{
  \item{x}{An instance of the rpa class (obtained as output from
	   RPA.pointestimate)} }

\value{An 'ExpressionSet' object.}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE. See
http://www.cis.hut.fi/projects/mi/software/RPA/}

\author{Leo Lahti <leo.lahti@iki.fi>}

\examples{

# Not run:

#require(RPA)
#require(affy)
#require(affydata)
#data(Dilution)

## Compute RPA for specific probesets
#sets <- geneNames(Dilution)[1:2]    
#rpa.results <- RPA.pointestimate(Dilution,sets)

## Coerce the rpa object into an ExpressionSet
#eset <- rpa2eset(rpa.results)

}

\keyword{ utilities }

