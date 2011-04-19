##setClass("TrioSet", contains="SnpCallSetPlus", representation(isBiparental="logical"))


setClass("LikSet",
	 contains="LogRatioSet",
	 representation(loglik="array",
			range.index="integer"),
	 prototype = prototype(
	 new("VersionedBiobase",
	     versions=c(classVersion("LogRatioSet"), LikSet="1.0.0"))))
setGeneric("loglik", function(object) standardGeneric("loglik"))
setGeneric("loglik<-", function(object, value) standardGeneric("loglik<-"))
setMethod("loglik", signature(object="LikSet"), function(object) object@loglik)
setReplaceMethod("loglik", signature(object="LikSet", value="array"), function(object,value){
	object@loglik <- value
	object
})
setReplaceMethod("loglik", signature(object="LikSet", value="numeric"), function(object,value){
	object@loglik <- value
	object
})
setGeneric("range.index", function(object) standardGeneric("range.index"))
setMethod("range.index", signature(object="LikSet"), function(object) fData(object)$range.index)
setMethod("[", "LikSet", function(x, i, j, ..., drop=FALSE){
	x <- callNextMethod(x, i, j, ..., drop=drop)
	loglik(x) <- loglik(x)[ , i, j, , ..., drop=drop]
	x
})

##setGeneric("isBiparental", function(object) standardGeneric("isBiparental"))
####setMethod("isBiparental", "TrioSet", function(object) object@isBiparental)
##setMethod("initialize", "TrioSet", function(.Object,
##					    isBiparental, ...){
##	.Object <- callNextMethod(.Object, ...)
##	.Object@isBiparental <- isBiparental
##})


