#' @importFrom parallel mclapply
#' @import affy
#' @importFrom BiocGenerics normalize

.onAttach <- function(lib, pkg)
{
   packageStartupMessage("\nRPA Copyright (C) 2008-2016 Leo Lahti.\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it under the FreeBSD open source license.\n")
}
