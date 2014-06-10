library(gridExtra)
library(ggplot2)
library(reshape)

# piclist, minsizes, accs.list, ntree, fs

k <- 2
fils <- c("Classification-comparisons.aff2.RData")

dfm <- NULL
pics <- list()
for (f in fils) {

  load(f)

  if (f == "Classification-comparisons-ENSG.RData") {
    maintext <- "Ensembl"
  } else if (f == "Classification-comparisons.aff2.RData") {
    maintext <- "Affymetrix"
  }

  accs <- accs.list[[as.character(k)]][, c("rma", "rpa.online", "rpa.online2", "frma", "frmab")]
  colnames(accs) <- c("RMA", "RPA", "RPA2", "fRMA", "fRMA.batch")

  df <- data.frame(accs); 
  df$fold <- 1:nrow(accs)
  df <- melt(df, id = c("fold"))
  df$method <- df$variable

  df$method <- factor(df$method, levels = c("RMA", "fRMA", "RPA", "RPA2", "fRMA.batch"))

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

  lims <- range(accs) 

  pics[[f]] <- ggplot(df) + aes(x = method, y = value) + theme_set(theme_bw(20)) + ylab("Classification accuracy") + xlab("") + opts(title = maintext, axis.text.x = theme_text(angle = 25, size = 20)) + opts(legend.position = "none") + geom_boxplot() + scale_y_continuous(limits = lims) +  opts(axis.text.x = theme_text(vjust = 0.5))

}

print(pics[[1]])

