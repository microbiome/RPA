\name{rpa.list-class}
%\Rdversion{1.1}
\docType{class}
\alias{rpa.list-class}
\alias{plot,rpa.list-method}

\title{Class "rpa.list"}
\description{Class for the RPA package.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rpa.list", ...)}.
}

\section{Slots}{
	 \describe{

An extended list. Contains the following information for one probeset.
 
    \item{\code{.Data}:}{Object of class \code{"list"}}

    \item{d }{A vector (probeset x arrays). Specifies the estimated
	      'true' underlying differential gene expression signal over
	      the arrays (vs. the reference array 'cind') for the
	      investigated probeset. Note that the reference array is not
	      included.}

    \item{mu.real }{Deviation between the signal shape and the real signal level. Here d is the sum of mu.real and shape; affinities are calculated around d = mu.real + shape. Note that in Lahti11 paper d denotes the shape; we will polish the notation later.}
	    
    \item{tau2 }{Contains a vector that gives the estimated variance
	           for each probe in the investigated probeset.}


  \item{affinity}{A list. Each element corresponds to a
		 probeset, and contains a vector that gives the
		 estimated affinity effect for each probe in that
		 probeset.}
		 
    \item{cind }{Specifies which of the arrays in abatch was used as the
	         reference for computing probe-level differential
	         expression.}

    \item{set }{Probeset name.}

    \item{data }{Preprocessed probe-level data on which the model was
		 fitted.}

  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "rpa.list")}: ... }
	 }
}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE 2011.}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\examples{showClass("rpa.list")}
\keyword{classes}
