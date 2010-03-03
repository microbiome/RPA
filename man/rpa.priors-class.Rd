\name{rpa.priors-class}
\Rdversion{1.1}
\docType{class}
\alias{rpa.priors-class}

\title{Class "rpa.priors"}

\description{Class for defining prior parameters in the RPA package.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rpa.priors", ...)}.

An extended list with elements alpha, beta, d.

\describe{

\item{alpha }{A list. Each element is a list that corresponds to one
		 probeset and specified the shape parameters of the
		 inverse Gamma distribution for the probe-specific
		 variance priors.}
						 
\item{beta }{A list. Each element is a list that corresponds to one
	     probeset and specified the scale parameters of the
	     inverse Gamma distribution for the probe-specific
	     variance priors.}
					     
\item{d }{Not implemented. Can be later used to set priors for d.}
}
}
%\section{Slots}{
%	 \describe{
%    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
%  }
%}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{No methods defined with class "rpa.priors" in the signature.}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE. See
http://www.cis.hut.fi/projects/mi/software/RPA/}

\author{Leo Lahti <leo.lahti@iki.fi>}

\examples{showClass("rpa.priors")}
\keyword{classes}
