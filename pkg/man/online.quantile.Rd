\name{online.quantile}
\alias{online.quantile}
\title{online.quantile
Quantile normalization tools for online preprocessing. Estimate quantiles for quantile normalization based on subset of the data (random, or specified by the user).}
\usage{
  online.quantile(abatch, n)
}
\arguments{
  \item{abatch}{AffyBatch}

  \item{n}{Numeric: number of random samples to use to
  define quantile basis. Vector: specify samples to be used
  in quantile basis calculation.}
}
\value{
  "online.quantile": AffyBatch; "qnorm.basis.online": a
  vector containing the basis for quantile normalization.
}
\description{
  online.quantile Quantile normalization tools for online
  preprocessing. Estimate quantiles for quantile
  normalization based on subset of the data (random, or
  specified by the user).
}
\details{
  "online.quantile": Ordinary quantile normalization is
  exhaustively memory-consuming in alrge data sets. Then
  the quantiles can be calculated based on subset of the
  data to allow efficient normalization. This function can
  also be used to investigate effect of subset size to
  convergence of the quantile
  estimates;"qnorm.basis.online": sweeps through the data
  in batches to calculate the basis for quantile
  normalization (average over sorted profiles).
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

