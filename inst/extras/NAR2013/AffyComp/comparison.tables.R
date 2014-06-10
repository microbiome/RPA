# this calculates mean ranking for the methods across the different tests
# test*.tab is manually copied from http://affycomp.jhsph.edu/

library(quantreg)

ranked.methods <- list()
ranked.list <- list()
# NOTE: test95.tab and test133.tab manually copy-pasted from 
# result site:
# http://affycomp.jhsph.edu/AFFY2/TABLES.hgu/0.html

# Comparison tests, names copypasted from affycomp
comp.names <- c("Median SD",
	        "null log-fc IQR",
		"null log-fc 99.9%",
		"Signal detect R2",
		"Signal detect slope",
		"low.slope",
		"med.slope",
		"high.slope",
		"Obs-intended-fc slope",
		"Obs-(low)int-fc slope", 
		"low AUC",
		"med AUC",
		"high AUC",
		"weighted avg AUC")


for (testfile in c("test95.tab.gz", "test133.tab.gz") ) {
  tab <- read.csv(testfile, sep = "\t", header = FALSE)

  # Pick method names
  mets <- gsub(" ", "", as.character(tab[,2]))

  # Polish result table
  res <- abs(apply(tab[, -c(1,2)], 1, as.numeric))
  colnames(res) <- mets
  rownames(res) <- comp.names
  res <- res[, !(duplicated(t(res)) & duplicated(colnames(res)))]

  # Calculate (abs.) deviation from ideal score for each method
  ideal.result <- res[,1]
  mat <- abs(apply(res, 2, function(x) {x - ideal.result}))

  # Calculate ranking for the methods
  ranked <- apply(mat, 1, function(x){match(1:ncol(mat), order(x))})
  rownames(ranked) <- colnames(mat)
  ranked.list[[testfile]] <- ranked

  # Results table
  resmat <- round(rbind(mat, rowMeans(ranked)),2)
  rownames(resmat) <- c(rownames(mat), "Average ranking")
  ranked.methods[[testfile]] <- resmat

  require(MASS)
  o <- order(resmat["Average ranking",])
  res.tab <- t(rbind(res[, o], resmat["Average ranking", o]))
  colnames(res.tab) <- c(colnames(res.tab)[1:14], "Average ranking")
  #write.matrix(res.tab, file = paste("results-", testfile, sep = ""))
  write.table(res.tab, file = paste("results-", testfile, sep = ""), quote = FALSE, col.names = TRUE, sep = " ")

  
  latex.table(res.tab, file=paste("Table-", testfile, sep = ""), digits = 2, lines.page = 40)

}

plist <- list()
comptabs <- list()
for (testfile in c("test95.tab.gz", "test133.tab.gz") ) {
  tab <- read.csv(testfile, sep = "\t", header = FALSE)

  # Pick method names
  mets <- gsub(" ", "", as.character(tab[,2]))

  # Polish result table
  res <- abs(apply(tab[, -c(1,2)], 1, as.numeric))
  colnames(res) <- mets
  rownames(res) <- comp.names
  #res <- res[, !(duplicated(t(res)) & duplicated(colnames(res)))]


  if (testfile == "test95.tab.gz") {chip <- "HG-U95Av2"}
  if (testfile == "test133.tab.gz") {chip <- "HG-U133A_tag"}
  mat2 <- res
  for (k in 1:3) {
      mat2[k,] <- 1 - mat2[k,]
  }

  mets <- c("RPA/leo.lahti", "RMA/rafa", "MAS_5.0/rafa")
  #mets <- c("RPA/leo.lahti", "rma/mboareto", "MAS_5.0/rafa")
  rrr <- mat2[, match(mets, colnames(mat2))]
  colnames(rrr) <- c("RPA", "RMA", "MAS5")
  df <- data.frame(rrr)
  comptabs[[testfile]] <- df
  rownames(df) <- 1:14
  #df$method <- rownames(df)
  library(reshape)
  dfm <- melt(t(df))
  names(dfm) <- c("Method", "Test", "value")
  dfm$Method <- factor(dfm$Method, levels = c("RPA", "RMA", "MAS5"))

  library(ggplot2)
  theme_set(theme_bw(20)); 
  p <- ggplot(dfm, aes(x = Test, y = value, group = Method, fill = Method)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_grey() + xlab("AffyComp Test") + ylab("Score") + scale_x_continuous(breaks=seq(from=1,to=14,by=1)) + ylim(limit = c(0, 1)) + ggtitle(chip)
  plist[[testfile]] <- p

}

library(gridExtra)
postscript("~/papers/rpa-online12/NAR/Lahti_Fig1.eps", onefile = FALSE, horizontal = FALSE, width = 14, height = 5, paper = "special", family = "Times")
grid.arrange(plist[[1]], plist[[2]], nrow = 1)
dev.off()



##################
#dset <- "test95.tab"; names(sort(ranked.methods[[dset]]["Average ranking",])[2:29])
#dset <- "test133.tab"; names(sort(ranked.methods[[dset]]["Average ranking",])[2:20])



#par(mfrow=c(2,1))
#for (dset in c("test95.tab", "test133.tab")) {
#  plot(sort(ranked.methods[[dset]]["Average ranking", ]), main = dset)
#  x <- match("RPA/leo.lahti", names(sort(ranked.methods[[dset]]["Average rankin#g",] )))
#  y <- sort(ranked.methods[[dset]]["Average ranking", "RPA/leo.lahti"])
#  points(x,y,pch=19)
#}


#> sort(ranked.methods[["test133.tab"]]["Average ranking", c("RPA/leo.lahti", "R#MA/rafa", "MAS_5.0/rafa")])
#RPA/leo.lahti      RMA/rafa  MAS_5.0/rafa 
#        33.71         35.64         54.64 
#> sort(ranked.methods[["test95.tab"]]["Average ranking", c("RPA/leo.lahti", "RM#A/rafa", "MAS_5.0/rafa")])
#RPA/leo.lahti      RMA/rafa  MAS_5.0/rafa 
#        32.71         33.57         44.71 

#tli.table <- xtable(res.tab)
#digits(tli.table)[c(2,6)] <- 0
#print(tli.table)
#print(tli.table, type="latex")

