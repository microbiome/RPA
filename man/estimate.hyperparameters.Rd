\name{estimate.hyperparameters}
\alias{estimate.hyperparameters}
\title{Hyperparameter estimation}
\description{Hyperparameter estimation}
\usage{
estimate.hyperparameters(sets = NULL, priors =  list(alpha = 2, beta = 1), batches, 
cdf = NULL, quantile.basis, bg.method = "rma", epsilon = 1e-2, cind = 1, 
load.batches = NULL, save.hyperparameter.batches = NULL, mc.cores = 1, 
verbose = TRUE, normalization.method = "quantiles")
}
\arguments{
  \item{sets}{Probesets to handle. All probesets by default.}
  \item{priors}{User-defined priors}
  \item{batches}{Data batches for online learning}
  \item{cdf}{CDF probeset definition file}
  \item{quantile.basis}{Basis for quantile normalization}
  \item{bg.method}{Background correction method}
  \item{epsilon}{Convergence parameter}
  \item{cind}{Reference array index}
  \item{load.batches}{Load preprocessed data with these file identifiers.}
  \item{save.hyperparameter.batches}{Save hyperparameters for each batch into files
    with this identifier.}
  \item{mc.cores}{Number of cores for parallel computation}
  \item{verbose}{Print progress information}
  \item{normalization.method}{Normalization method}
}
\details{Sweeps through the batches. Updates hyperparameters at each batch.}
\value{
  \item{alpha }{Hyperparameter alpha (same for all probesets)}
  \item{betas }{Hyperparameter beta (probe-specific)}
  \item{variances }{Probe-specific variances (beta/alpha)}

}
\references{See citation("RPA").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\seealso{}
\examples{
#
}
\keyword{ utilities }
