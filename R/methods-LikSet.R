setMethod("range.index", signature(object="LikSet"), function(object) fData(object)$range.index)
setMethod("[", "LikSet", function(x, i, j, ..., drop=FALSE){
	x <- callNextMethod(x, i, j, ..., drop=drop)
	loglik(x) <- loglik(x)[ , i, j, , ..., drop=drop]
	x
})

setMethod("loglik", signature(object="LikSet"), function(object) object@loglik)
setReplaceMethod("loglik", signature(object="LikSet", value="array"), function(object,value){
	object@loglik <- value
	object
})
setReplaceMethod("loglik", signature(object="LikSet", value="numeric"), function(object,value){
	object@loglik <- value
	object
})
