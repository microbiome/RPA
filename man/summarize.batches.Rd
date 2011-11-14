\name{summarize.batches}
\alias{summarize.batches}
\title{Summarize batches}
\description{Probeset summarization in batches}
\usage{
summarize.batches(sets, variances, batches, load.batches, mc.cores = 1,
cdf = NULL, bg.method = "rma", normalization.method = "quantiles", verbose = FALSE, quantile.basis)
}
\arguments{
  \item{sets}{Probesets to summarize}
  \item{variances}{Precalculated probe-specific variances}
  \item{batches}{Data batches for online learning}
  \item{load.batches}{Batch identifier string for precalculated batch data}
  \item{mc.cores}{Number of cores for parallel computation}
  \item{cdf}{CDF for alternative probeset definitions}
  \item{bg.method}{Background correction method}
  \item{normalization.method}{Normalization method}
  \item{verbose}{Print progress information}
  \item{quantile.basis}{Basis for quantile normalization}
}
\details{Sweeps through the batches. Summarizes the probesets within
  each batch based on the precalculated model parameters using point
  estimates of the model.}
\value{An ExpressionSet object}
\references{See citation("RPA").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\seealso{}
\examples{
#
}
\keyword{ utilities }

