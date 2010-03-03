#setMethod(f = "plot", signature("rpa.list"),
#        function(x, y){
#        }
#)


plot.rpa.list <- function (x, y, ...) {

      set <- x$set
        d <- x$d
       s2 <- x$sigma2		  

      par(mfrow = c(2,1))
      barplot(s2, main = paste(set,"/ Probe-specific variances (sigma2)"),
              ylab = "Variance", xlab = "Probes", names.arg = 1:length(s2), las = 2)              
      barplot(d, main = paste(set,"/ Differential expression estimate (d)"),
              ylab = "Dif. exp. signal", xlab = "Arrays", las = 2)	      

}

#plot.rpa <- function (x, y, set = NULL...) {
#  # y: probeset
#  # x: rpa.object to plot  
#
#  # if set not specified, plot the first one
#  set <- ifelse(is.null(set), x$sets[[1]], set)
#  rpa.plot(set, x, ...)
#}


