\name{rpa.complete}
\alias{rpa.complete}
\title{rpa.complete
RPA preprocessing, also returns probe parameters}
\usage{
  rpa.complete(abatch = NULL, sets = NULL, epsilon = 0.01,
    tau2.method = "robust", d.method = "fast",
    verbose = FALSE, bg.method = "rma",
    normalization.method = "quantiles.robust", cdf = NULL,
    cel.files = NULL, cel.path = NULL,
    probe.parameters = list(), mc.cores = 1,
    summarize.with.affinities = FALSE)
}
\arguments{
  \item{abatch}{An AffyBatch object.}

  \item{sets}{Probesets for which RPA will be computed.}

  \item{epsilon}{Convergence tolerance. The iteration is
  deemed converged when the change in all parameters is <
  epsilon.}

  \item{tau2.method}{Optimization method for tau2
  (probe-specific variances). This parameter is denoted by
  tau^2 in the vignette and manuscript

  "robust": (default) update tau2 by posterior mean,
  regularized by informative priors that are identical for
  all probes (user-specified by setting scalar values for
  alpha, beta). This regularizes the solution, and avoids
  overfitting where a single probe obtains infinite
  reliability. This is a potential problem in the other
  tau2 update methods with non-informative variance priors.
  The default values alpha = 2; beta = 1 are used if alpha
  and beta are not specified.

  "mode": update tau2 with posterior mean

  "mean": update tau2 with posterior mean

  "var": update tau2 with variance around d. Applies the
  fact that tau2 cost function converges to variance with
  large sample sizes.}

  \item{d.method}{Method to optimize d.

  "fast": (default) weighted mean over the probes, weighted
  by probe variances The solution converges to this with
  large sample size.

  "basic": optimization scheme to find a mode used in Lahti
  et al. TCBB/IEEE; relatively slow; this is the preferred
  method with small sample sizes.}

  \item{verbose}{Print progress information during
  computation.}

  \item{bg.method}{Specify background correction method.
  Default: "rma". See bgcorrect.methods() for other
  options.}

  \item{normalization.method}{Specify quantile
  normalization method. Default: "pmonly". See
  normalize.methods(Dilution) for other options.}

  \item{cdf}{Specify an alternative CDF environment.
  Default: none.}

  \item{cel.files}{List of CEL files to preprocess.}

  \item{cel.path}{Path to CEL file directory.}

  \item{probe.parameters}{A list, each element
  corresponding to a probe set. Each probeset element has
  the following optional elements: affinity (affinity),
  tau2 (variance), alpha (shape prior), betas (scale
  prior). Each of these elements contains a vector over the
  probeset probes, specifying the probe parameters
  according to the RPA model. If variance is given, it
  overrides the priors. Can be also used to set
  user-specified priors for the model parameters. Not used
  tau2.method = "var". The prior parameters alpha and beta
  are prior parameters for inverse Gamma distribution of
  probe-specific variances. Noninformative prior is
  obtained with alpha, beta -> 0.  Not used with
  tau2.method 'var'. Scalar alpha and beta specify an
  identical inverse Gamma prior for all probes, which
  regularizes the solution. Can be also specified as lists,
  each element corresponding to one probeset. Can also
  include quantile.basis}

  \item{mc.cores}{Number of cores for parallelized
  processing.}

  \item{summarize.with.affinities}{Use affinity estimates
  in probe summarization step. Default: FALSE.}
}
\value{
  List with preprocessed expression matrix, corresponding
  probe parameters, AffyBatch and CDF
}
\description{
  rpa.complete RPA preprocessing, also returns probe
  parameters
}
\details{
  RPA preprocessing function. Gives an estimate of the
  probeset-level mean parameter d of the RPA model, and
  returns these in an expressionSet object. The choices
  tau2.method = "robust" and d.method = "fast" are
  recommended. With small sample size and informative
  prior, d.method = "basic" may be preferable. For very
  large expression data collections, see rpa.online
  function.
}
\examples{
# eset <- rpa(abatch)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("RPA")
}
\keyword{methods}

