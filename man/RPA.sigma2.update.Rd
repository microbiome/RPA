\name{RPA.sigma2.update}

\Rdversion{1.1}

\alias{RPA.sigma2.update}

\title{Updating probe-specific variances.}

\description{Updates the probe-specific variance parameters sigma2 given
data S and estimated mean d, assuming the prior parameters alphahat,
beta}

\usage{RPA.sigma2.update(d, S = S, alphahat, beta, sigma2.method="var")}

\arguments{

  \item{d }{A vector. Estimated 'true' signal underlying the noisy
  probe-level observations.}
  
  \item{S }{Matrix of probe-level observations for a single probeset:
  samples x probes.}

  \item{alphahat }{T/2 + alpha, where T = nrow(S) (sample size), and
  alpha is the shape parameter for inverse Gamma prior.}

  \item{beta }{Scale parameter for inverse gamma prior.}
  
  \item{sigma2.method }{
      Optimization method for sigma2 (probe-specific variances).
      
      "basic": optimization using user-specified alpha, beta priors.
	       Default: alpha, beta = 1e-6.
		 
      "var": utilizes the fact that the cost function converges to
      	     variance with large sample sizes. Default method.	       
   }
}

\value{A vector of probe-specific variances.}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE, to appear. See
http://www.cis.hut.fi/projects/mi/software/RPA/ }

\author{Leo Lahti <leo.lahti@tkk.fi>}

\seealso{RPA.sigma2.update.variance, RPA.sigma2.update.fast}

\keyword{ iteration }
\keyword{ internal }
