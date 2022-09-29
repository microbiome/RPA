#' @import affy
#' @importFrom BiocGenerics normalize
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom parallel mclapply
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom stats dchisq
#' @importFrom stats dgamma
#' @importFrom stats rgamma
#' @importFrom stats optim
#' @importFrom stats na.omit
#' @importFrom utils read.csv
#' @importFrom utils sessionInfo

.onAttach <- function(lib, pkg)
{
   packageStartupMessage("\nRPA Copyright (C) 2008-2021 Leo Lahti. See http://microbiome.github.io/\n")
}
