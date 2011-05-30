\name{sample.probeset}
\alias{sample.probeset}
\title{Toydata generator for probeset data.}
\description{Generate random probeset with varying probe-specific affinities and variances.}
\usage{
sample.probeset(P = 10, n = 20, shape = 1, scale = 1, mu.real = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{Number of probes.}
  \item{n}{Number of samples.}
  \item{shape, scale}{Shape and scale parameters of the inverse Gamma function used to generate the probe-specific variances.}
  \item{mu.real}{Absolute signal level of the probeset.}
}
\details{The toy data generator follows distributional assumptions of the RPA model and allows quantitative estimation of model accuracy with different options, noise levels and sample sizes. Probeset-level summary estimate is obtained as mu.real + d.}
\value{
 A list with the following elements:
 \item{dat }{Probeset data: probes x samples}
 \item{variance }{Probe variances.}
 \item{affinity }{Probe affinities.}
 \item{d }{Probeset-level signal shape.}
 \item{mu.real }{Probeset-level absolute signal level.}

}
\references{See citation("RPA").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\examples{
real <- sample.probeset(P = 10, n = 20, shape = 1, scale = 1, mu.real = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }

