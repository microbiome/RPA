\name{rpa-class}
\Rdversion{1.1}
\docType{class}
\alias{rpa-class}
\alias{[,rpa-method}
\alias{[[,rpa-method}
\alias{show,rpa-method}

\title{Class "rpa"}
\description{Class for the RPA package.}

\section{Objects from the Class}{
Returned by \code{\link{RPA.pointestimate}} function.
Objects can be created by calls of the form \code{new("rpa", ...)}.
}

\section{Slots}{
    \describe{

   Contains the following information for analyzed probesets and data:
  
   \item{\code{.Data}:}{Object of class \code{"list"}}

   \item{d }{A matrix of probesets x arrays. Specifies the estimated
  	'true' underlying gene expression signal over the arrays for
  	each investigated probeset.}

   \item{mu.real }{Numeric. Used for technical purposes. Indicates the
   difference between original signal levels and the signal levels in
   the differential signal domain used to learn the probe-specific
   variance parameters. See model description in vignette for details.}
      
  \item{tau2 }{A list. Each element corresponds to a probeset, and
		 contains a vector that gives the estimated variance
		 (noise level) for each probe in that probeset.}
    
  \item{affinity}{A list. Each element corresponds to a
		 probeset, and contains a vector that gives the
		 estimated affinity effect for each probe in that
		 probeset.}
				  
   \item{cind }{Specifies which of the arrays in abatch was used as a
		reference for computing probe-level differential
		expression.}
	
   \item{sets }{A character vector listing the investigated
                probesets.}

   \item{cdf }{Alternative CDF that was used in the analysis.}

   \item{data }{Preprocessed probe-level data on which the model was
     fitted.}

   \item{abatch }{The associated affybatch object.}
  }

}

\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "rpa")}: ... }
    \item{[[}{\code{signature(x = "rpa")}: ... }
    \item{show}{\code{signature(x = "rpa")}: ... }
  }
}


\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE. http://www.roihu.info/publications/preprints/TCBB10.pdf}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\examples{showClass("rpa")}
\keyword{classes}
