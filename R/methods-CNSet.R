setMethod("inheritance", "CNSet", function(object){
	stopifnot("inheritance" %in% fvarLabels(object))
	fData(object)$inheritance
})
setMethod("mendelian", "CNSet", function(object){
	stopifnot("mendelian" %in% fvarLabels(object))
	fData(object)$mendelian
})
setMethod("informative", "CNSet", function(object){
	stopifnot("informative" %in% fvarLabels(object))
	fData(object)$mendelian
})

setMethod("axis.limits", "CNSet", function(object){
	if(nrow(object) > 1) {
		warning("multiple snps in object.  using the first.")
		object <- object[1,]
	}
	x <- c(log2(A(object)), log2(B(object)))
	axis.limit <- range(x, na.rm=TRUE) + c(-1,1)*0.5
})

setMethod("gtConfidence", "CNSet", function(object){
	1-exp(-callsConfidence(object)/1000)
})
