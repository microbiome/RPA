# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. 
# All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' rpa2eset
#' Coerce 'rpa' object into an 'ExpressionSet'
#'
#' @param x An instance of the rpa class (obtained as output from RPA.pointestimate)
#'
#' @details An instance of 'rpa' class contains differential gene expression estimates in the variable 'd'. The function 'rpa2eset' coerces this into an ExpressionSet object to allow downstream analysis of the results using standard R/BioC tools for gene expression data.
#'
#' @return An 'ExpressionSet' object.
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords utilities

rpa2eset <- function(x) {

  # x: an rpa object                                        
  # Coerces expression values in the rpa object into an ExpressionSet object

  require(affy)
  new("ExpressionSet", 
	assayData = list(exprs = x$d), 
	phenoData = phenoData(x), 
	annotation = annotation(x), 
        protocolData = protocolData(x), 
	experimentData = experimentData(x))

}

