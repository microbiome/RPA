library(RPA)

library(ALLMLL)
data(MLL.A)
t.rma <- system.time(eset.rma <- rma(MLL.A))
t.rpa <- system.time(eset.rpa <- rpa(MLL.A))
cors <- c()
for (i in rownames(eset.rma)) {cors[[i]] <- cor(eset.rma[i, ], eset.rpa[i, ])}

cors <- c()
for (i in rownames(eset.rma)) {cors[[i]] <- cor(eset.rma[i, ], eset.rpa[i, ])}

#pdf("~/tmp/tmp.pdf")
hist(cors, 100)
#dev.off()

rpa.results <- RPA.pointestimate(MLL.A)

pdf("~/tmp/tmp.pdf")
par(mfrow=c(3,3))
for (k in 1:9) {
  set <- names(sort(cors))[[k]] #rownames(eset.rpa)[[k]]
  dat <- plot.rpa(rpa.results, set = set, external.signal = eset.rma[set,], main = set)
}
dev.off()


rpa.results$d[1,]
eset.rpa[1,]
