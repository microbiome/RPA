\name{initialize.priors}
\Rdversion{1.1}
\alias{initialize.priors}

\title{Initialize prior parameter object for RPA.}

\description{Creates a template for prior parameters using
	     non-informative (or user-specified) priors. After
	     creating the template, probe-specific priors can be set
	     by modifying the default values. Current version supports
	     priors are only for the probe-specific variances
	     (i.e. priors alpha, beta).}

\usage{initialize.priors(abatch, sets, alpha = 1e-6, beta = 1e-6, d = NULL)}

\arguments{

\item{abatch}{An AffyBatch object.}
	    
\item{sets}{Vector listing the probesets for which the prior template
	    is created.}

\item{alpha}{Default template values for the alpha prior}

\item{beta}{Default template values for the beta prior}

\item{d}{Not used, included for later compatibility.}

}

\details{Probesets that do not have predefined priors are analyzed
	 using noninformative alpha, beta (1e-6). See the source code
	 for RPA.iteration for details.}

\value{An instance of 'rpa.priors' class.}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE, to appear. See
http://www.cis.hut.fi/projects/mi/software/RPA/}

\author{Leo Lahti <leo.lahti@tkk.fi>}

\examples{

require(affy)
require(affydata)
data(Dilution)

# List probesets to investigate
sets = geneNames(Dilution)[1:5]             

# Create a template for prior parameters using default priors
my.priors = initialize.priors(Dilution, sets, alpha = 1e-6, beta = 1e-6)

# Modify the template to provide user-specified prior for one of the
# probes
# high value implies an unreliable probe
set = sets[[1]]
probe.idx = 5
my.priors[[set]]$beta[[probe.idx]] = 10  
my.priors[[set]]$alpha[[probe.idx]] = 10

# Run RPA using the predefined priors
# NOTE: priors are only used with sigma2.method = "basic"
rpa.results <- RPA.pointestimate(Dilution,sets,
	       			 priors=my.priors,
	       			 sigma2.method = "basic", 
				 d.method = "basic")


# This toy example shows the probe reliability values for the probeset
# where the user-specified prior was set for one of the probes
barplot(rpa.results$sigma2[[set]])

}

\keyword{ methods }

