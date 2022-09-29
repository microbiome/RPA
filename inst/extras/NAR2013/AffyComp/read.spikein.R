library(affycomp)
filename = "frma-rwa.csv" 
cdfName = "hgu133a"
remove.xhyb=TRUE

fs <- list.files("~/Rpackages/RPA/AffyComp/affycomp/R/", full.names = TRUE)
for (f in fs) {source(f)}

#    cdfName <- match.arg(cdfName)
    s <- read.csv(filename, check.names = FALSE, row.names = 1)
    samplenames <- colnames(s)
    samplenames <- sub("\\.gz$", "", samplenames, ignore.case = TRUE)
    samplenames <- sub("\\.Z$", "", samplenames, ignore.case = TRUE)
    samplenames <- sub("\\.cel$", "", samplenames, ignore.case = TRUE)
    colnames(s) <- samplenames
    if (cdfName == "hgu95a") {
        spikein.phenodata <- getData("spikein.phenodata")
        pd <- spikein.phenodata
    }
    if (cdfName == "hgu133a") {
        hgu133a.spikein.phenodata <- getData("hgu133a.spikein.phenodata")
        pd <- hgu133a.spikein.phenodata
    }
    s <- s[, sampleNames(pd)]
    s <- new("ExpressionSet", exprs = as.matrix(s), phenoData = pd)
    s <- exprset.log(s)
    if (remove.xhyb & cdfName == "hgu133a") 
        s <- remove.hgu133a.xhyb(s)
    return(s)
