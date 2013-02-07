\name{rpa.fit-class}
\Rdversion{1.1}
\docType{class}
\alias{rpa.fit-class}
%\alias{[,rpa.fit-method}
%\alias{[[,rpa.fit-method}
%\alias{show,rpa.fit-method}

\title{Class "rpa.fit"}
\description{Class for the RPA package.}

\section{Objects from the Class}{
Returned by \code{\link{rpa.fit}} function.
Objects can be created by calls of the form \code{new("rpa.fit", ...)}.
}

\section{Slots}{
    \describe{

   Contains the following information for analyzed probesets and data:
  
   \item{\code{.Data}:}{Object of class \code{"list"}}

   \item{mu }{A matrix of probes x samples. Gives the summary estimate
     corresponding to the 'true' underlying signal over the samples.}

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
  \item{data }{Original probe-level data.}
  \item{alpha, beta }{Inverse Gamma prior parameters for variance (tau2).}


  }
}

\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}

%\section{Methods}{
%  \describe{
%    \item{[}{\code{signature(x = "rpa.fit")}: ... }
%    \item{[[}{\code{signature(x = "rpa.fit")}: ... }
%    \item{show}{\code{signature(x = "rpa.fit")}: ... }
%  }
%}


\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE 2011. http://www.roihu.info/publications/preprints/TCBB10.pdf}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\examples{showClass("rpa.fit")}
\keyword{classes}
