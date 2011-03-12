
# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

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

