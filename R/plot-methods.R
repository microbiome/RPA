# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


plot.rpa.list <- function (x, y, ...) {

      set <- x$set # probe set
       mu <- x$d # 'Absolute' expression level
       sd <- sqrt(x$tau2) # probe standard deviation	 

      par(mfrow = c(2,1))
      barplot(sd, main = paste(set,"/ Probe noise (standard deviation)"),
              ylab = "Standard deviation", xlab = "Probes", names.arg = 1:length(sd), las = 2)

      barplot(mu, main = paste(set,"/ expression "),
              ylab = "Expression signal", xlab = "Arrays", las = 1) 

}


