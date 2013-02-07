# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' rpa.plot
#' Plot RPA results and probe-level data for a specified probeset.
#'
#' @param dat Original data: probes x samples.
#' @param rpa.fit.object An instance of the 'rpa.fit' class.
#' @param toydata.object Optional. Output from sample.probeset toydata generator function. Can be used to compare (toy)data with known ground truth to RPA estimates from rpa.fit.object.
#' @param highlight.probes Optionally highlight some of the probes (with dashed line)
#' @param pcol Color for probe signal visualization.
#' @param mucol Color for summary estimate.
#' @param ecol Color for external signal.
#' @param cex.lab Label size adjustment parameters.
#' @param cex.axis Axis size adjustment parameters.
#' @param cex.main Title size adjustment parameters.
#' @param cex.names Names size adjustment parameters.
#' @param external.signal Plot external signal on the probeset. For instance, an alternative summary estimate from another preprocessing methods
#' @param main Title text.
#' @param plots "all": plot data and summary, noise and affinity; "data": plot data and summary
#' @param ... Other parameters to pass for plot function.
#'
#' @details Plots the preprocessed probe-level observations, estimated probeset-level signal, and probe-specific variances. It is also possible to highlight individual probes and external summary measures.
#'
#'@return Used for its side-effects. Returns probes x samples matrix of probe-level data plotted on the image.
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords methods

rpa.plot <- function (dat, rpa.fit.object = NULL, toydata.object = NULL, highlight.probes = NULL, pcol = "darkgrey", mucol = "black", ecol = "red", cex.lab = 1.5, cex.axis = 1, cex.main = 1, cex.names = 1, external.signal = NULL, main = "", plots = "all", ...) {

  # rpa.fit.object = estimated; toydata.object = real; plots = "toydata"; highlight.probes = NULL; pcol = "darkgrey"; mucol = "black"; ecol = "red"; cex.lab = 1.5; cex.axis = 1; cex.main = 1; cex.names = 1; external.signal = NULL; main = ""

  # If no model given, calculate fit a new model on the data
  if (is.null(rpa.fit.object)) { rpa.fit.object <- rpa.fit(dat) }

  # Override given data with the one in rpa.fit object
  dat <- rpa.fit.object$data
  
  # number of probes
  Np <- nrow(dat)

  mu <- rpa.fit.object$mu
  sd <- sqrt(rpa.fit.object$tau2)
  af <- rpa.fit.object$affinity
  
  # image limits
  ylims <- range(c(as.vector(dat), mu))

  if (plots == "all") { par(mfrow = c(3, 1)) }
  if (plots == "toydata") { par(mfrow = c(2, 2)) }
  
  # expression figure
  plot(c(1,2,3), type = 'n',
       xlim = c(1, ncol(dat)),
       ylim = ylims,
       xlab = "Samples",
       ylab = "Signal",
       cex.lab = cex.lab,
       cex.axis = cex.axis,
       cex.main = cex.main,
       main = main,
       las = 1, xaxt = 'n', ...)
  #axis(1, at = 1:ncol(dat), labels = colnames(dat), las = 2, cex.axis = 0.5)
  #axis(1, at = 1:ncol(dat), las = 2, cex.axis = 0.5)

  for (i in 1:Np) { lines(dat[i, ], col = pcol, lwd = 2) }

  if (!is.null(highlight.probes)) {
    lines(dat[highlight.probes, ], lty = 2, lwd = 2) 
  }

  # Plot the summary
  lines(mu, col = mucol, lwd = 2)

  # Plot external signal
  if (!is.null(external.signal)) {
    lines(external.signal, col = ecol, lwd = 2)
  }

  if (plots == "all") {
  
    # Plot probe variances
    barplot(sd, main = "Stochastic probe effect (noise)",
          ylab = "Noise (tau)", xlab = "Probe index",
          names.arg = 1:length(sd),
          las = 1,
          cex.lab = cex.lab,
          cex.axis = cex.axis, cex.names = cex.names, cex.main = cex.main)

    # Plot probe affinities
    barplot(af, main = "Fixed probe effect (affinity)",
          ylab = "Affinity (mu)", xlab = "Probe index",
          names.arg = 1:length(af),
          las = 1,
          cex.lab = cex.lab, 
          cex.axis = cex.axis, cex.names = cex.names, cex.main = cex.main)  

  } else if (plots == "toydata") {

    estimated <- rpa.fit.object
    real <- toydata.object

    plot(real$d + real$mu.real, estimated$mu, 
    	          main = "Signal", xlab = "Real", ylab = "Estimated")
    abline(0,1)
    barplot(rbind(real$affinity, estimated$affinity), 
                  beside = TRUE, main = "Probe affinity", xlab = "Probe index", ylab = "Affinity")

    tab <- rbind(real = real$variance, estimated = estimated$tau2)
    barplot(tab, beside = TRUE, 
            main = "Probe variance", 
            xlab = "Probe index", ylab = "Variance", legend = TRUE)

  }
  
  rpa.fit.object

}
