RPA.preprocess <-
function (abatch, cind) {

	#print("Preprocess raw data")
	#Get raw expression matrix and perform quantile normalization
	q<-log2(exprs(normalize.AffyBatch.quantiles(abatch, type="separate")))
	colnames(q)<-colnames(exprs(abatch))

	#Select control sample randomly if not given
	if (!cind) {cind<-sample(ncol(q),1)}

	#Compute fold changes against the control sample
	fcmat <- q[,-cind]-q[,cind]
	colnames(fcmat) <- colnames(q)[-cind]

	return(list(fcmat = fcmat, cind = cind))
}

