\name{online.quantile}
\alias{online.quantile}
\alias{qnorm.basis.online}

\title{Quantile normalization tools for online preprocessing.}
\description{Estimate quantiles for quantile normalization based on subset of the data (random, or specified by the user).}
\usage{
online.quantile(abatch, n)
qnorm.basis.online(cel.files, bg.method = "rma", cdf = NULL, save.batches = NULL, batch.size = 2)
}
\arguments{
  \item{abatch}{AffyBatch}
  \item{n}{Numeric: number of random samples to use to define quantile basis. Vector: specify samples to be used in quantile basis calculation.}
  \item{cel.files}{List of CEL files to process}
  \item{bg.method}{Background correction method.}
  \item{cdf}{Optional. CDF file for alternative probeset definitions.}
  \item{save.batches}{Optional. Save results for each batch into file
  with identifier 'save.batches-1.RData' etc.}
  \item{batch.size}{Batch size.}
}
\details{"online.quantile": Ordinary quantile normalization is
  exhaustively memory-consuming in alrge data sets. Then the quantiles
  can be calculated based on subset of the data to allow efficient
  normalization. This function can also be used to investigate effect of
  subset size to convergence of the quantile estimates;
  "qnorm.basis.online": sweeps through the data in batches to calculate
  the basis for quantile normalization (average over sorted profiles).}
\value{"online.quantile": AffyBatch; "qnorm.basis.online": a vector
  containing the basis for quantile normalization.}
\references{See citation("RPA").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\seealso{RPA.preprocess}
\examples{
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }



