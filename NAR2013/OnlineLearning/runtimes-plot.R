load("runtimes.RData")

rt <- sapply(runtimes, function (x) {x[["elapsed"]]})

library(ggplot2)
theme_set(theme_bw(20))
p <- qplot(ns, rt/(60^2))
p <- p + xlab("Sample size (arrays)") + ylab("Time (h)") + geom_line() #+ opts(axis.text.x = theme_text(size = 20))

postscript("runtimes.eps", onefile = TRUE)
print(p)
dev.off()

pdf("runtimes.pdf")
print(p)
dev.off()
