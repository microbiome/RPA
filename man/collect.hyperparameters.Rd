\name{collect.hyperparameters}
\alias{collect.hyperparameters}
\title{collect.hyperparameters}
\usage{
  collect.hyperparameters(batches, unique.run.identifier, save.batches.dir, batch.size)
}
\arguments{
  \item{batches}{batch list}
  \item{unique.run.identifier}{unique run identifier}
  \item{save.batches.dir}{Batch data file directory}
  \item{batch.size}{Batch size}
}
\value{
  List; one element for each parameter. With updated parameter values at each batch.
}
\description{
  Use this function to investigation parameter evolution during online-learning.
}
\examples{
# hpe <- collect.hyperparameters(batches, unique.run.identifier, save.batches.dir, batch.size)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("RPA")
}
\keyword{utilities}

