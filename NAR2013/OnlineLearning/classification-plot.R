library(gridExtra)
library(ggplot2)
library(reshape)

# piclist, minsizes, accs.list, ntree, fs

for (k in c(2)) {

#k <- 2 # 2 : minsize 5
fils <- c("Classification-comparisons.RData", 
          "Classification-comparisons-ENSG.RData")

dfm <- NULL
pics <- list()
lims <- c(0.82, 0.90)
theme_set(theme_bw(20))
for (f in fils) {

  load(f)

  if (f == "Classification-comparisons-ENSG.RData") {
    maintext <- "Ensembl"
  } else if (f == "Classification-comparisons.RData") {
    maintext <- "Affymetrix"
  }

  accs <- accs.list[[as.character(k)]][, c("mas", "rma", "rpa.online", "frma", "frmab")]
  colnames(accs) <- c("MAS5", "RMA", "RPA", "fRMA", "fRMA.batch")
  #o <- order(apply(accs, 2, mean))
  #accs <- accs[,o]
  df <- data.frame(accs); 
  df$fold <- 1:nrow(accs)
  df <- melt(df, id = c("fold"))
  #df <- melt(df)
  df$method <- df$variable

  df$method <- factor(df$method, levels = c("MAS5", "RMA", "fRMA", "RPA", "fRMA.batch"))

  #df$label <- gsub("Online.RPA", "Online-RPA", as.character(df$method))
  #df$label <- gsub("fRMA.batch", "fRMA-batch", as.character(df$label))

  # Testing if row method is less accurate than col method
  pvals <- matrix(NA, nrow = ncol(accs), ncol = ncol(accs))
  rownames(pvals) <- colnames(pvals) <- colnames(accs)
  for (method1 in colnames(accs)) {
    for (method2 in colnames(accs)) {
      pvals[method1, method2] <- wilcox.test(accs[, method1], accs[, method2], "less", paired = TRUE)$p.value
    } 
  }    
  print(round(pvals, 5))

  write.table(round(pvals,5), sep = "\t", file = paste("Classification-comparisons.", maintext, k, ".tab", sep = ""))

  pics[[f]] <- ggplot(df) + aes(x = method, y = value) + theme_set(theme_bw(20)) + ylab("Classification accuracy") + xlab("") + opts(title = maintext, axis.text.x = theme_text(angle = 25, size = 20)) + opts(legend.position = "none") + geom_boxplot() + scale_y_continuous(limits = lims) +  opts(axis.text.x = theme_text(vjust = 0.5))

}

minsizes <- c(1,2,3)
#postscript(paste("Classification-comparisons.", minsizes[[k]], ".eps", sep = ""), width = 14, height = 7, onefile = FALSE, horizontal = FALSE, paper = "special", family = "Times")
postscript(paste("~/papers/rpa-online12/NAR/Lahti_Fig3.eps", sep = ""), width = 14, height = 7, onefile = FALSE, horizontal = FALSE, paper = "special", family = "Times")
grid.arrange(pics[[1]], pics[[2]], nrow = 1)
dev.off()

pdf(paste("~/papers/rpa-online12/NAR/Lahti_Fig3.pdf", sep = ""), width = 14, height = 7, family = "Times")
grid.arrange(pics[[1]], pics[[2]], nrow = 1)
dev.off()


}



#ggplot(df) + aes(x = method, y = value) + geom_line()

