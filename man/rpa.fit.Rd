\name{rpa.fit}
\alias{rpa.fit}
\title{rpa.fit}
\description{Fits the RPA model, including estimation of probe-specific
  affinity parameters.}
\usage{
rpa.fit(dat, cind = 1, epsilon = 1e-2, alpha = NULL, beta = NULL,
tau2.method = "robust", d.method = "fast")
}

\arguments{
  \item{dat}{Original data: probes x samples.}
  \item{cind}{Index of reference array.}
  \item{epsilon }{Convergence tolerance. The iteration is deemed
  converged when the change in all parameters is < epsilon.}
  
  \item{alpha, beta }{Priors for inverse Gamma distribution of
  probe-specific variances. Noninformative prior is obtained with alpha,
  beta -> 0.  Not used with tau2.method 'var'. Scalar alpha and beta
  are specify equal inverse Gamma prior for all probes to regularize the
  solution. The defaults depend on the method.}

  \item{tau2.method }{

    Optimization method for tau2 (probe-specific variances).

	"robust": (default) update tau2 by posterior mean,
		regularized by informative priors that are identical
		for all probes (user-specified by
		setting scalar values for alpha, beta). This
		regularizes the solution, and avoids overfitting where
		a single probe obtains infinite reliability. This is a
	        potential problem in the other tau2 update
	        methods with non-informative variance priors. The
		default values alpha = 2; beta = 1 are
	        used if alpha and beta are not specified.
   
        "mode": update tau2 with posterior mean

	"mean": update tau2 with posterior mean
	
	"var": update tau2 with variance around d. Applies the fact
               that tau2 cost function converges to variance with
               large sample sizes. 

		}

  \item{d.method }{
    Method to optimize d.


        "fast": (default) weighted mean over the probes, weighted by
		probe variances The solution converges to this with
		large sample size.

        "basic": optimization scheme to find a mode used in Lahti et
        	 al. TCBB/IEEE; relatively slow; preferred 
		 with small sample size.
                
	      }

}
\details{First learns a point estimate for the RPA model in terms of
  differential expression values w.r.t. reference sample. After this,
  probe affinities are estimated by comparing original data and
  differential expression shape, and setting prior assumptions
  concerning probe affinities.}
\value{

  \item{mu}{Fitted signal in original data: mu.real + d}
  \item{mu.real}{Shifting parameter of the reference sample}
  \item{tau2}{Probe-specific stochastic noise}
  \item{affinity}{Probe-specific affinities}
  \item{data}{Probeset data matrix}
  \item{alpha, beta}{prior parameters}

}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE 2011. See
http://www.roihu.info/publications/preprints/TCBB10.pdf}

\author{Leo Lahti \email{leo.lahti@iki.fi}}


\note{Affinity estimation is not part of the original RPA procedure in
TCBB/IEEE 2011 paper. It is added here since estimates of the absolute
levels are often needed in microarray applications. Note that affinity
parameters are unidentifiable in the model if no prior assumptions are
given. We assume that affinity effects are zero on average, but allow
some flexibility through probe-specific weights.}

\seealso{rpa, RPA.pointestimate, estimate.affinities}
\examples{
## res <- rpa.fit(dat, cind, epsilon, alpha, beta, tau2.method, d.method, affinity.method)
}
\keyword{ utilities }

