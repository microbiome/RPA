# Expression matrix for publication
#write.table(round(emat, 3), file = "HGU133A-Lukk2010-Affy-RPA-ExpressionMatrix-20130216.tab", sep = "\t", quote = FALSE, row.names = FALSE)
#probe.parameters.table <- probetable(res$probe.parameters)
#probe.parameters.list  <- probe.parameters.tolist(probe.parameters.table)
#probe.parameters.table <- probetable(probe.parameters.list)

# Hyperparameter table for publication
#probe.parameters.rounded[, c("alpha", "betas", "tau2", "affinity")] <- round(probe.parameters.table[, c("alpha", "betas", "tau2", "affinity")], 3)
#write.table(probe.parameters.rounded, file = "HGU133A-Lukk2010-Affy-RPA-ProbeParams-20130216.tab", sep = "\t", quote = FALSE, row.names = FALSE)

#load("~/data/RPA/Lukk-RPA.RData")
#emat3 <- emat.rpa[rownames(emat),  colnames(emat)]
#rs <- sample(prod(dim(emat)), 1e5)
#plot(emat[rs], emat3[rs], pch = '.')

# load("LukkOutput20130216.RData")
#require(RPA2);
#probe.parameters <- probetable(probe.parameters)
#save(emat, probe.parameters, file = "Online-RPA-Lukk2010-Affy.RData", compress = "xz")

probetable <- function (hyper.parameters) {

  # Arrange probe-level parameters into a table
  sets <- unlist(sapply(1:length(hyper.parameters$affinities), function (k) {rep(names(hyper.parameters$affinities)[[k]], length(hyper.parameters$affinities[[k]]))}))
  inds <- unlist(sapply(1:length(hyper.parameters$affinities), function (k) {1:length(hyper.parameters$affinities[[k]])}))
  df <- data.frame(list(probeset = sets, probe.index = inds))
  df$alpha <- rep(hyper.parameters$alpha, nrow(df)) # probe hyperparameter alpha_j (same for all probes; function of sample size T)
  df$beta <- unlist(hyper.parameters$betas)     	  # probe hyperparameter beta_j
  df$tau2 <- unlist(hyper.parameters$variances) 	  # probe variance tau2_j (calculated from alpha, beta)
  df$mu <- unlist(hyper.parameters$affinities)  	  # probe affinity mu_j   (has a prior from tau2)
  df

}
get.probe.parameters <- function (unique.run.identifier, save.batches.dir = ".", mode = "list") {

  # Final estimated hyperparameters (alpha, betas, variances)
  load(paste(save.batches.dir, "/", unique.run.identifier, "-RPA-hyperparameters.RData", sep = "")) # hyper.parameters
  load(paste(save.batches.dir, "/", unique.run.identifier, "-RPA-affinities.RData", sep = "")) # affinities
  hyper.parameters$affinities <- affinities

  hp <- hyper.parameters
  if (mode == "table") {
    hp <- probetable(hyper.parameters)
  } 

  hp 

}


#length(unlist(pmindex(abatch)))

#mask.inds <- xy2indices(xy[,1], xy[,2], abatch = abatch)
#mask.sets <- c()
#for (set in featureNames(abatch)) {
#  if (any(unlist(pmindex(abatch, set)) %in% mask.inds)) {
#    print(set)
#    mask.sets <- c(mask.sets, set)
#  }
#}

#########################################
#indices2xy(unlist(pmindex(abatch, set))[[1]], abatch = abatch)
                                        #xy2indices(xy[,1], xy[,2], abatch = abatch)
# pmindex(abatch, set)

#print(dim(xy))                                        
#print(pmindex(abatch, set))

#eset <- rma(abatch)
#set <- report$removedprobe$probeset[[1]]
#barplot(rbind(exprs(eset.filtered)[set,], exprs(eset)[set,]), beside = TRUE)

#xy <- getsnpprobe(chip = abatch@cdfName)
#xy <- getsnpprobe(chip = "HGU133Plus2")
#xy<- getsnpprobe(chip='HGU133Plus2', criteria='het > 0.5')
#cdf <- 
#if (!is.null(cdf)) { 

# loc > 0 and loc < 10
# het > 0
# rsid in (1,3,4,5)
# loc   Location of allele on the probe, upper value is probe's length,
#       lower limit could be negative if the allele is longer than probe
# het	Heterozygosity value
# rsid	Snp id




fs <- list.files("../RPA/R", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}

set.inds <- get.set.inds(batches[[1]][1:2])

abatch <- ReadAffy(filenames = cel.files)
abatch.bg <- bg.correct(abatch, bg.method, destructive = TRUE)
ab.q2 <- do.call(affy:::normalize, c(alist(abatch.bg, "quantiles"), normalize.param = list()))
qmat2 <- log2(pm(ab.q2))

sets <- names(set.inds)
set <- sets[[1]]
s2 <- variances[[set]]
dat <- qmat2[set.inds[[set]],]

set.summary <- colSums(dat/s2)/sum(1/s2)
cor(eset.rma[set,s], set.summary[s])
cor(emat[set,s], set.summary[s])

rpa.online <- function (
 	        cel.path = NULL,
               cel.files = NULL, 
                    sets = NULL,
                  priors = list(alpha = 2, beta = 1),
                 epsilon = 1e-2, 
                    cind = 1,
                 verbose = FALSE,
               bg.method = "rma",
	       normalization.method = "quantiles",
                     cdf = NULL, 
              batch.size = 10, 
	  quantile.basis = NULL) 
{

  }

  if (is.null(cel.files) && !is.null(cel.path)) {
    cel.files <- list.celfiles(cel.path, full.names = TRUE)
  }

  #################################################################
  
  # NOTE: list CEL files in random order to avoid biases!
  # FIXME: add this
  cel.files.sampled <- cel.files # sample(cel.files)

  # Get probe position indices
  abatch <- ReadAffy(filenames = cel.files[1:2], compress=getOption("BioC")$affy$compress.cel) 
  # Set alternative CDF environment if given
  if (!is.null(cdf)) { abatch@cdfName <- cdf }    

  # Check names and number for the investigated probesets
  # if my.sets not specified, take all sets in abatch
  if (is.null(sets)) { sets <- geneNames(abatch) } 
  Nsets <- length( sets ) 

  # Retrieve probe positions
  # The indices correspond to the output of pm()
  pN <- probeNames(abatch)
  set.inds <- split(1:length(pN), pN)[sets] # pNList
  
  ###############################################################

  if (normalization.method == "quantiles") {
    # Online-estimation of the basis for quantile normalization
    message("Calculating the basis for quantile normalization")
    quantile.basis <- qnorm.basis.online(cel.files.sampled, bg.method, batch.size, cdf)
  }
  
  ###############################################################
  
  # Split CEL file list into batches
  # TODO: list CEL files in random order to avoid biases!
  # TODO: Define batching already before quantile.basis calculation?
  batches <- get.batches(cel.files.sampled, batch.size)
 
  ###############################################################

  hyps <- estimate.hyperparameters(priors, sets, set.inds, batches, cdf, quantile.basis, bg.method, normalization.method, epsilon)
  alpha <- hyps$alpha  
  betas <- hyps$betas

  ###############################################################

  # Get final estimated variances for each probeset based on hyperparameter posteriors
  variances <- lapply(betas, function (beta) {beta/alpha})
  names(variances) <- names(betas) 
  
  # Now hyperparameters for probe variances have been estimated
  # Get probeset-level signal estimate by online sweep with the fixed variances
  emat <- array(NA, dim = c(length(sets), length(cel.files.sampled)))
  rownames(emat) <- sets
  colnames(emat) <- cel.files.sampled # sapply(strsplit(cel.files, "/"), function (x) {x[[length(x)]]})

  for (i in 1:length(batches)) {

    message(paste("Summarizing batch", i, "/", length(batches)))

    # CEL file IDs for this batch
    batch.cels <- batches[[i]]
    
    # Get background corrected, quantile normalized, and logged probe-level matrix
    # Do NOT calculate probe-level diff.exp here any more, only needed in variance estimation.
    q <- get.probe.matrix(cels = batch.cels, cdf, quantile.basis, bg.method, normalization.method)

    # Get probes x samples matrices for each probeset
    q <- lapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) })
    names(q) <- sets

    for (set in sets) {
      # print(which(set == sets)/length(sets))    
      # Get summary estimate using the posterior variance
      emat[set, batch.cels] <- d.update.fast(q[[set]], variances[[set]])
    }
  }

  # Return the arrays into original order
  emat <- emat[, cel.files]
  
  # Remove path from CEL names for compactness
  # colnames(emat) <- sapply(strsplit(colnames(emat), "/"), function (x) { x[[length(x)]] })

  # Coerce expression values in the rpa object into an ExpressionSet object
  # and return expression set
  new("ExpressionSet", assayData = list(exprs = emat)) 
  
}


                                        # normalize(cel.files)
#abatch <- ReadAffy(filenames = cel.files)
#abatch.bg2 <- abatch <- do.call(affy:::bg.correct,
#                  c(alist(abatch, method = "rma"), bgcorrect.param = list()))



#system.time(apply(pm(abatch.bg), 2, function(x) {o <- rank(x); log2(rowMeans(apply(pm(abatch.bg), 2, sort)))[o]}))
#system.time(log2(pm(do.call(affy:::normalize, c(alist(abatch.bg, "quantiles"), normalize.param = list())))))
                                        #rownames(qmat) <- rownames(pm(abatch.bg))

#cor(order(qmat2[,1]), qmat[,1])


#rownames(qmat2)
#rownames(qmat)

#pmb1 <- pm(abatch.bg)[,1]


#cors <- diag(cor(qmat, qmat2[, colna

#set.quantiles <- function (mat, quantile.basis) {
#  # replace smallest value with the smallest in the given quantile.basis etc.
#  apply(mat, 2, function(x) { 
#    x[order(x)] <- sort(quantile.basis);
#    x
#  })
#}

                                        # quantile.basis <- quantile.basis.online(cel.files, bg.method, batch.size, cdf)mes(qmat)]))
#print(cors)
#o1 <- order(pm(abatch.bg)[,1])
#   eset <- computeExprSet(afbatch, summary.method = summary.method, 
#        pmcorrect.method = pmcorrect.method, ids = summary.subset, 
#        summary.param = summary.param, pmcorrect.param = pmcorrect.param)
#emat.rpa.online <- exprs(res)

########################################################

# STANDARD RPA
#emat.rpa.standard <- exprs(rpa(cel.files = cels, priors = list(alpha = 2, beta = 1), verbose = TRUE))
# RMA
#emat.rma <- exprs(rma(ReadAffy(filenames = cels)))
#rsample <- sample(prod(dim(emat.rma)), 1e3)
#tab <- cbind(RMA = emat.rma[rsample], RPA.online = emat.rpa.online[rsample], RPA.standard = emat.rpa.standard[rsample])
#print(cor(tab, use = "pairwise.complete.obs"))


#(r2.colsums - r.colsums^2 / ((nrow(R) + 1)))  *(1/nrow(R))
#r.colsums <- colSums(R)
#r2.colsums <- colSums(R^2)
#beta + (r2.colsums - r.colsums^2 / ((nrow(R) + 1)))/2
#N <- nrow(R)
#print("1-")
#print((N/2)*colSums(R^2)/(N-1))
#print("2-")
#print((N/2)*apply(R,2,var))
#print("3-")
#print(colSums(R^2)/2)
#print("3-")






update.beta.EM <- function (dat, alpha, beta, th) {

  # FIXME: combine with update.beta

  # Estimate s2 mode with EM-type procedure, initialize with prebvious s2
  s2 <- s2.update(dat, alpha, beta, s2.init = beta/alpha, th = th) # FIXME move conv. param. to arguments

  # Return beta; calculated from the mode of s2, which is s2mode <- beta/alpha
  alpha * s2 # <- beta

}


# Generate and fit toydata, learn hyperparameters
P <- 11   # number of probes
N <- 1000 # number of arrays
real <- sample.probeset(P = P, n = N, shape = 3, scale = 1, mu.real = 4)
dat <- real$dat # probes x samples

res <- rpa.online.test(dat, priors = list(alpha = 1, beta = 1), batch.size = 500)
res$alpha
res$beta
est.var <- res$beta/res$alpha
res$d

plot(sqrt(real$variance), sqrt(est.var), main = cor(sqrt(real$variance), sqrt(est.var)))




set.seed(11122)
nsets <- 50

source("known.abatch.R")
