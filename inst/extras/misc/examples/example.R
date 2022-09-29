library(RPA)
#if (!requireNamespace("BiocManager", quietly=TRUE))
    #install.packages("BiocManager")
#BiocManager::install("ALLMLL")
library(ALLMLL)
data(MLL.A)
abatch <- MLL.A

#######################

rpa.results <- RPA.pointestimate(abatch)
#k <- k+1
#set <- rpa.results$sets[[k]]
set <- "200046_at"
#plot(rpa.results, set = set, plots = "all", cex.axis = 2)
postscript("example.eps")
par(mar=c(4, 4.5, 2, 1) + 0.1)
plot(rpa.results, set = set, plots = "all", cex.axis = 1.5, cex.main = 2, cex.names = 2, cex.lab = 2)
dev.off()

############################

#library(ggplot2)
#dat <- as.data.frame(rpa.results$data[pmindices,])
#dat[["probe"]] <- as.character(1:nrow(dat))

#datm <- melt(dat)
#colnames(dat) <- c("probe", "array", "signal")
#mu <- data.frame(list(array = names(x$d[1,]), signal = x$d[1,]))

#p <- ggplot(datm, aes(x = array, y = signal, group = probe)) + geom_line(colour = "gray")
#p + ggplot(mu) + geom_line(aes(x = array, y = signal), colour = "black") 
print(p)


