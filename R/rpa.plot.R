#' @title rpa.plot
#' @description Plot RPA results and probe-level data for a specified probeset.
#'
#' @param x Output from rpa.complete function
#' @param set probeset
#' @param highlight.probes mark probes for highlight
#' @param pcol probe color
#' @param mucol probeset signal color
#' @param ecol external signal color
#' @param external.signal external signal to be plotted on top 
#' @param main title
#' @param plots plot type
#' @param ... other arguments to be passed
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

rpa.plot <- function (x, set, highlight.probes = NULL, pcol = "darkgrey", mucol = "black", ecol = "red", external.signal = NULL, main = NULL, plots = "all", ...) {
                                      
  # get the associated affybatch
  # Use alternative CDF environment if given                                        

  abatch <- x$abatch
  if (!is.null(x$cdf)) { abatch@cdfName <- x$cdf }

  # Find probe (pm) indices for this set
  pmindices <- which(probeNames(abatch) == set)

  # Pick all model parameters
      tau2 <- x$probe.parameters$tau2[[set]]
  affinity <- x$probe.parameters$affinity[[set]]
        mu <- exprs(x$eset)[set,]
       dat <- x$probedata[pmindices,]
 
  if (is.null(main)) {
    main <- paste("Probe signals and the summary estimate (", set, ")", sep = "")
  }

  rpaplot(dat, mu, tau2, affinity, main = main, plots = plots, external.signal = external.signal, pcol = pcol, ecol = ecol, mucol = mucol, ...)

  dat
}






#' rpaplot
#' Plot RPA results and probe-level data for a specified probeset.
#'
#' @param dat Background-corrected and normalized data: probes x samples.
#' @param mu probeset signal
#' @param tau2 probe variances
#' @param affinity probe affinities
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

rpaplot <- function (dat, mu = NULL, tau2 = NULL, affinity = NULL, highlight.probes = NULL, pcol = "darkgrey", mucol = "black", ecol = "red", cex.lab = 1.5, cex.axis = 1, cex.main = 1, cex.names = 1, external.signal = NULL, main = "", plots = "all", ...) {

  sd <- sqrt(tau2)

  if (plots == "all") { par(mfrow = c(3, 1)) }
  if (plots == "toydata") { par(mfrow = c(2, 2)) }
  
  # Plot the probes
  tmp <- probeplot(dat, highlight.probes = NULL, pcol = "darkgrey", mucol = "black", ecol = "red", cex.lab = 1.5, cex.axis = 1, cex.main = 1, cex.names = 1, external.signal = NULL, main = "", ...)

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
    barplot(affinity, main = "Fixed probe effect (affinity)",
          ylab = "Affinity", xlab = "Probe index",
          names.arg = 1:length(affinity),
          las = 1,
          cex.lab = cex.lab, 
          cex.axis = cex.axis, cex.names = cex.names, cex.main = cex.main)  

  } 
  
  NULL

}







#' probeplot
#' Plot RPA results and probe-level data for a specified probeset.
#'
#' @param dat Background-corrected and normalized data: probes x samples.
#' @param highlight.probes Optionally highlight some of the probes (with dashed line)
#' @param pcol Color for probe signal visualization.
#' @param hcol Color for probe highlight
#' @param cex.lab Label size adjustment parameters.
#' @param cex.axis Axis size adjustment parameters.
#' @param cex.main Title size adjustment parameters.
#' @param cex.names Names size adjustment parameters.
#' @param main Title text.
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

probeplot <- function (dat, highlight.probes = NULL, pcol = "darkgrey", hcol = "red", cex.lab = 1.5, cex.axis = 1, cex.main = 1, cex.names = 1, main = "", ...) {

  # highlight.probes = NULL; pcol = "darkgrey"; hcol = "red"; cex.lab = 1.5; cex.axis = 1; cex.main = 1; cex.names = 1; main = "";

  # image limits
  ylims <- range(c(as.vector(dat)))

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

  for (i in 1:nrow(dat)) { 
    lines(dat[i, ], col = pcol, lwd = 2) 
  }

  for (probe in highlight.probes) {
    lines(dat[probe, ], lwd = 1, col = hcol) 
  }

  NULL

}


