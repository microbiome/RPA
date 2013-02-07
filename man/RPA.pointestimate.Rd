\name{RPA.pointestimate}
\alias{RPA.pointestimate}
\title{RPA.pointestimate}
\usage{
  RPA.pointestimate(abatch, sets = NULL, myseed = 101,
    priors = NULL, epsilon = 0.01, cind = 1,
    tau2.method = "robust", d.method = "fast",
    verbose = TRUE, bg.method = "rma",
    normalization.method = "quantiles.robust", cdf = NULL)
}
\arguments{
  \item{abatch}{An AffyBatch object.}

  \item{sets}{Specifies the probesets for which RPA
  estimates will be computed. Default: all probe sets.}

  \item{myseed}{Specifies the random seed.}

  \item{priors}{Optional list containing hyperparameters
  alpha and beta of the inverse Gamma prior of the
  probe-specific variances; alpha is a scalar, common for
  all probes and probeset; beta is a list where each
  element is a vector corresponding to one probeset,
  specifying beta for each probe. Can be used to set
  user-specified priors for the model parameters. Not
  applicable for tau2.method = "var". Noninformative prior
  is obtained with alpha, beta approaches 0. Not used with
  tau2.method 'var'. Can be used to regularize the solution
  with small sample size, or in batch-wise online-updates
  with large sample size. If priors are not provided for
  certain probesets (NULL), default priors are used.}

  \item{epsilon}{Convergence tolerance. The iteration is
  deemed converged when the change in the d parameter is
  smaller than epsilon.}

  \item{cind}{Specifies which array in abatch is used as a
  reference in computing probe-level differential
  expression.}

  \item{tau2.method}{Optimization method for tau2
  (probe-specific variances). "robust": (default) update
  tau2 by posterior mean, regularized by informative priors
  that are identical for all probes (user-specified by
  setting scalar values for alpha, beta). This regularizes
  the solution and avoids overfitting where a single probe
  obtains infinite reliability. This is a potential problem
  in the other tau2 update methods with non-informative
  variance priors. The default values alpha = 2; beta = 1
  are used if alpha and beta are not specified; "mode":
  update tau2 with posterior mean; "mean": update tau2 with
  posterior mean; "var": update tau2 with variance around
  d. Applies the fact that tau2 cost function converges to
  variance with large sample size}

  \item{d.method}{Method to optimize d: "fast": (default)
  weighted mean over the probes, weighted by probe
  variances The solution converges to this with large
  sample size; "basic": optimization scheme to find a mode
  used in Lahti et al. TCBB/IEEE; relatively slow; this is
  the preferred method with small sample sizes.}

  \item{verbose}{Print progress information during
  computation.}

  \item{bg.method}{Specify background correction method.
  Default: "rma". See bgcorrect.methods() for other
  options.}

  \item{normalization.method}{Specify quantile
  normalization method. Default: "pmonly". See
  normalize.methods(Dilution) for other options.}

  \item{cdf}{Specify an alternative CDF environment.}
}
\value{
  An instance of class 'rpa'. This is an extended list
  containing the following elements: d: A matrix of
  probesets x arrays. Specifies the estimated 'true'
  underlying differential gene expression signal over the
  arrays (vs. the reference array 'cind') for each
  investigated probeset. Note that the reference array is
  not included; tau2: A list. Each element corresponds to a
  probeset, and contains a vector that gives the estimated
  variance for each probe in that probeset. This
  corresponds to the parameter tau^2 in the vignette and
  manuscript; cind Specifies which of the arrays in the
  abatch (the affybatch object to be analyzed) has been
  used as the reference for computing probe-level
  differential expression; affinity: Probe affinity
  effects; sets: A character vector listing the
  investigated probesets.
}
\description{
  Computing point estimate for the model parameters for all
  probe sets
}
\details{
  Calculates RPA estimates of probe reliability and
  differential expression between the user-specified
  reference array (cind) and the other arrays in the data
  set. The model assumes P observations for each transcript
  target (i.e. a probeset) with Gaussian noise which is
  specific for each probe (variance is specified by tau2).
  The mean (affinity) parameters of the Gaussian noise
  model cancel out in calculating probe-level differential
  expression.  RPA.pointestimate gives a point estimate for
  d and tau2. The 'prior' parameter is not applicable with
  tau2.method = "var". The d.method = "fast" is recommended
  with large sample size. tau2.method = "robust" and
  d.method = "fast" are recommended. With small sample size
  and informative priors, d.method = "basic" may be
  preferable, with large sample size d.method = "fast"
  should have considerable speedups and comparable
  accuracy.
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
\seealso{
  rpa.plot, rpa, set.priors, rpa2eset, RPA.preprocess,
  AffyBatch, rpa.fit, estimate.affinities
}
\keyword{methods}
\keyword{utilities}

