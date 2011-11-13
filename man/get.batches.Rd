\name{get.batches}
\alias{get.batches}
\title{Split data into batches}
\description{Splits data into batches}
\usage{
get.batches(items, batch.size, shuffle = FALSE)
}
\arguments{
  \item{items}{A vector of items to be splitted into batches.}
  \item{batch.size}{Batch size. The last batch may contain less elements
  than the other batches which have batch.size elements each.}
  \item{shuffle}{Split the elements randomly in the batches.}
}
%\details{}
\value{A list. Each element corresponds to one batch and contains a
  vector listing the elements in that batch.}
\references{See citation("RPA").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\seealso{RPA.preprocess}
\examples{
#
}
\keyword{ utilities }

