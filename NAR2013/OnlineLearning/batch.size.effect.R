# https://github.com/microbiome/microbiome/wiki/RPA

require(RPA)
set.seed(23)
cels <- sample(list.celfiles("CEL", full.names = T))
cels <- sample(cels, 300)

fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names = T); for (f in fs) {source(f)};

times <- list()
esets <- list()
bss <- c(20, 50, 100, 300)
for (bs in bss) {
  print(paste("Batch size: ", bs))
  times[[as.character(bs)]] <- system.time(eset <- rpa.online(cel.files = cels, save.batches.dir = "~/tmp/", batch.size = bs, mc.cores = 4))
  esets[[as.character(bs)]] <- eset
}


bs <- as.character(bss[[1]])
rs <- sample(prod(dim(exprs(esets[[bs]]$expressionSet))), 1e3)

mat <- NULL
for (k in 1:length(bss)) {
  mat <- cbind(mat, exprs(esets[[k]]$expressionSet)[rs])
}
colnames(mat) <- as.character(bss)
pairs(mat)
