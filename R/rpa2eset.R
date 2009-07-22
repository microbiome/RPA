rpa2eset <-
function(x) {
	# x: an rpa object
	# Coerces expression values in the rpa object into an ExpressionSet object
	require(affy)
	new("ExpressionSet", assayData = list(exprs=x$d))
}

