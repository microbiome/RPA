\name{RPA.iteration}
\alias{RPA.iteration}
\title{RPA.iteration
Estimating model parameters d and tau2.}
\usage{
  RPA.iteration(S, epsilon = 0.001, alpha = NULL,
    beta = NULL, tau2.method = "fast", d.method = "fast",
    maxloop = 1e+06)
}
\arguments{
  \item{S}{Matrix of probe-level observations for a single
  probeset: samples x probes.}

  \item{epsilon}{Convergence tolerance. The iteration is
  deemed converged when the change in all parameters is <
  epsilon.}

  \item{alpha}{alpha prior for inverse Gamma distribution
  of probe-specific variances. Noninformative prior is
  obtained with alpha, beta -> 0.  Not used with
  tau2.method 'var'. Scalar alpha and beta are specify
  equal inverse Gamma prior for all probes to regularize
  the solution. The defaults depend on the method.}

  \item{beta}{beta prior for inverse Gamma distribution of
  probe-specific variances. Noninformative prior is
  obtained with alpha, beta -> 0.  Not used with
  tau2.method 'var'. Scalar alpha and beta are specify
  equal inverse Gamma prior for all probes to regularize
  the solution. The defaults depend on the method.}

  \item{tau2.method}{Optimization method for tau2
  (probe-specific variances).

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

  \item{d.method}{Method to optimize d.  "fast": (default)
  weighted mean over the probes, weighted by probe
  variances The solution converges to this with large
  sample size.

  "basic": optimization scheme to find a mode used in Lahti
  et al. TCBB/IEEE; relatively slow; this is the preferred
  method with small sample sizes.}

  \item{maxloop}{Maximum number of iterations in the
  estimation process.}
}
\value{
  A list with the following elements: d: A vector.
  Estimated 'true' signal underlying the noisy probe-level
  observations.; tau2: A vector. Estimated variances for
  each measurement (or probe).
}
\description{
  RPA.iteration Estimating model parameters d and tau2.
}
\details{
  Finds point estimates of the model parameters d
  (estimated true signal underlying probe-level
  observations), and tau2 (probe-specific variances).
  Assuming data set S with P observations of signal d with
  Gaussian noise that is specific for each observation
  (specified by a vector tau2 of length P), this method
  gives a point estimate of d and tau2. Probe-level
  variance priors alpha, beta can be used with tau2.methods
  'robust', 'mode', and 'mean'.  The d.method = "fast" is
  the recommended method for point computing point
  estimates with large samples size.
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

