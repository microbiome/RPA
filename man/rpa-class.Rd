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
  	'true' underlying differential gene expression signal over the
  	arrays (vs. the control array 'cind') for each investigated
  	probeset. Note that the control array is not included.}
       
  \item{sigma2 }{A list. Each element corresponds to a probeset, and
		 contains a vector that gives the estimated variance
		 for each probe in that probeset.}
    
   \item{cind }{Specifies which of the arrays in abatch was used as a
		control in computing probe-level differential
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
al., TCBB/IEEE. See http://www.cis.hut.fi/projects/mi/software/RPA/}

\author{Leo Lahti <leo.lahti@iki.fi>}

\examples{showClass("rpa")}
\keyword{classes}
