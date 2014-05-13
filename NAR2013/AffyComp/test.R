#library(affycomp) ##load the package
#s <- read.spikein("frma-rwa.csv", cdfName = "hgu133a")

#ass <- assessSpikeIn(s,method.name="fRMA")

method.name="fRMA"
verbose = TRUE
#affycomp is a wrapper for assessAll and affycompTable. The re
    if (verbose) 
        cat("Performing 6 assessments that will take a few minutes")
#    tmp1 <- assessFC(s, method.name = method.name)


exprset <- s


    e <- exprs(exprset)
    pdata <- pData(exprset)

    if (ncol(e) == 59) {
        WHICHSPIKEIN <- "HGU95A"
    } else {
        if (ncol(e) == 42) 
            WHICHSPIKEIN <- "HGU133A"
        else stop("Not the right number of columns in expression matrix\n")
    }

    genenames <- colnames(pdata)
    N <- length(genenames)
    M <- nrow(e) - N
    intended <- array(0, dim = c(570, N, 2))
    observed <- matrix(0, 570, N)
    fc2 <- matrix(0, 570, 2)
    probs <- c(0, 25/M, 100/M, 0.25, 0.75, 1 - 100/M, 1 - 25/M, 
        1)
    quantiles <- matrix(0, 570, length(probs))
    pdata <- as.matrix(pdata)
    spikeindex <- match(genenames, rownames(e))
    roc <- vector(mode = "numeric", length = nrow(e))
    Count <- 0
    if (WHICHSPIKEIN == "HGU95A") {
        J <- 20
    } else {J <- 14 }


    for (i in 1:(J - 1)) {
        for (j in (i + 1):J) {

            i1 <- pdata[i, ]
            i2 <- pdata[j, ]

            if (!all(i1 - i2 == 0)) {
                Count <- Count + 1
                intended[Count, , 1] <- i1
                intended[Count, , 2] <- i2
                m <- e[, j] - e[, i]
print("JEI")
        quantiles[Count, ] <- quantile(m[-spikeindex], prob = probs)
                  


                fc2[Count, 1] <- sum(abs(m[-spikeindex]) >= 1)
                fc2[Count, 2] <- sum(abs(m[spikeindex]) >= 1)
                observed[Count, ] <- m[genenames]
                m <- sort(-abs(m))
                y <- rep(0, length(m))
                y[match(genenames, names(m))] <- 1
                roc <- roc + cumsum(y)
            }
        }
    }


    for (i in (J + 1):(2 * J - 1)) {
        for (j in (i + 1):(2 * J)) {
            i1 <- pdata[i, ]
            i2 <- pdata[j, ]
            if (!all(i1 - i2 == 0)) {
                Count <- Count + 1
                intended[Count, , 1] <- i1
                intended[Count, , 2] <- i2
                m <- e[, j] - e[, i]
                quantiles[Count, ] <- quantile(m[-spikeindex], 
                  prob = probs)
                fc2[Count, 1] <- sum(abs(m[-spikeindex]) >= 1)
                fc2[Count, 2] <- sum(abs(m[spikeindex]) >= 1)
                observed[Count, ] <- m[genenames]
                m <- sort(-abs(m))
                y <- rep(0, length(m))
                y[match(genenames, names(m))] <- 1
                roc <- roc + cumsum(y)
            }
        }
    }
    if (WHICHSPIKEIN == "HGU95A") {
        J1 <- 41
        J2 <- 58
    } else {
        J1 <- 29
        J2 <- 42
    }
    for (i in J1:(J2 - 1)) {
        for (j in (i + 1):J2) {
            i1 <- pdata[i, ]
            i2 <- pdata[j, ]
            if (!all(i1 - i2 == 0) & i != 54 & j != 54) {
                Count <- Count + 1
                intended[Count, , 1] <- i1
                intended[Count, , 2] <- i2
                m <- e[, j] - e[, i]
                quantiles[Count, ] <- quantile(m[-spikeindex], 
                  prob = probs)
                fc2[Count, 1] <- sum(abs(m[-spikeindex]) >= 1)
                fc2[Count, 2] <- sum(abs(m[spikeindex]) >= 1)
                observed[Count, ] <- m[genenames]
                m <- sort(-abs(m))
                y <- rep(0, length(m))
                y[match(genenames, names(m))] <- 1
                roc <- roc + cumsum(y)
            }
        }
    }
    intended <- intended[1:Count, , ]
    observed <- observed[1:Count, ]
    fc2 <- fc2[1:Count, ]
    quantiles <- quantiles[1:Count, ]
    quantiles <- colMeans(quantiles)
    tp <- roc/Count
    fp <- seq(along = tp) - tp
    colnames(observed) <- genenames
    dimnames(intended) <- list(NULL, genenames, NULL)
    names(quantiles) <- c("lowest", "lowest25", "lowest100", 
        "25", "75", "highest100", "highest25", "highest")
    intended.log.ratios <- log2(intended[, , 2]/intended[, , 
        1])
    x <- as.vector(intended.log.ratios)
    y <- as.vector(observed)
    Index <- as.vector(intended[, , 2]) <= 2 & as.vector(intended[, 
        , 1]) <= 2 & as.vector(intended[, , 2]) > 0 & as.vector(intended[, 
        , 1]) > 0
    N <- ncol(pdata)
    list(signal = intended, intended.log.ratios = intended.log.ratios, 
        observed.log.ratios = observed, quantiles = quantiles, 
        fc2 = fc2, tp = tp, fp = fp, area = c(a10 = mean(tp[fp < 
            10]/N), a15 = mean(tp[fp < 15]/N), a25 = mean(tp[fp < 
            25]/N), a100 = mean(tp[fp < 100]/N)), slope = lm(y ~ 
            x, subset = abs(x) < Inf)$coef[2], low.signal.slope = lm(y ~ 
            x, subset = Index)$coef[2], index.low.signal = Index, 
        what = "FC", method.name = method.name)
