\name{online.quantile}
\alias{online.quantile}
\title{Quantile normalization for online preprocessing.}
\description{Estimate quantiles for quantile normalization based on subset of the data (random, or specified by the user).}
\usage{
online.quantile(abatch, n)
}
\arguments{
  \item{abatch}{AffyBatch}
  \item{n}{Numeric: number of random samples to use to define quantile basis. Vector: specify samples to be used in quantile basis calculation.}
}
\details{Ordinary quantile normalization is exhaustively memory-consuming in alrge data sets. Then the quantiles can be calculated based on subset of the data to allow efficient normalization. This function can also be used to investigate effect of subset size to convergence of the quantile estimates.}
\value{AffyBatch}
\references{See citation("RPA").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\seealso{RPA.preprocess}
\examples{
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }

