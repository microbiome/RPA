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


