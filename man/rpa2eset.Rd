\name{rpa2eset}
\alias{rpa2eset}
\title{rpa2eset
Coerce 'rpa' object into an 'ExpressionSet'}
\usage{
  rpa2eset(x)
}
\arguments{
  \item{x}{An instance of the rpa class (obtained as output
  from RPA.pointestimate)}
}
\value{
  An 'ExpressionSet' object.
}
\description{
  rpa2eset Coerce 'rpa' object into an 'ExpressionSet'
}
\details{
  An instance of 'rpa' class contains differential gene
  expression estimates in the variable 'd'. The function
  'rpa2eset' coerces this into an ExpressionSet object to
  allow downstream analysis of the results using standard
  R/BioC tools for gene expression data.
}
\examples{
#
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("RPA")
}
\keyword{utilities}

