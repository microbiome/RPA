\name{summarize.batches}
\alias{summarize.batches}
\title{Summarize batches}
\description{Probeset summarization in batches}
\usage{
summarize.batches(sets = NULL, variances, batches, load.batches = FALSE,
mc.cores = 1, cdf = NULL, bg.method = "rma", normalization.method =
"quantiles", verbose = TRUE, quantile.basis, save.batches.dir = ".", unique.run.identifier = NULL, save.batches = FALSE, affinities)
}
\arguments{
  \item{sets}{Probesets to summarize}
  \item{variances}{Precalculated probe-specific variances}
  \item{batches}{Data batches for online learning}
  \item{load.batches}{Logical. Load precalculated data for the batches.}
  \item{mc.cores}{Number of cores for parallel computation}
  \item{cdf}{CDF for alternative probeset definitions}
  \item{bg.method}{Background correction method}
  \item{normalization.method}{Normalization method}
  \item{verbose}{Print progress information}
  \item{quantile.basis}{Basis for quantile normalization}
  \item{save.batches.dir}{Specify the output directory for temporary batch saves.}
  \item{unique.run.identifier}{Define identifier for this run for naming the temporary batch files. By default, a random id is generated.}
  \item{save.batches}{Save batches?}
  \item{affinities}{probe-specific affinities}

}
\details{Sweeps through the batches. Summarizes the probesets within
  each batch based on the precalculated model parameter point estimates.}
\value{Expression matrix: probesets x samples.}
\references{See citation("RPA").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\seealso{}
\examples{
#
}
\keyword{ utilities }

