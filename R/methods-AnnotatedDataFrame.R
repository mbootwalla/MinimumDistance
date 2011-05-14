annotatedDataFrameFromArray <- function(object, byrow=FALSE, ...){
	object <- object[, , 1, drop=TRUE]
	res <- Biobase:::annotatedDataFrameFromMatrix(object, byrow=byrow, ...)
	return(res)
}
setMethod("annotatedDataFrameFrom", signature(object="ff_array"),
	  annotatedDataFrameFromArray)
setMethod("annotatedDataFrameFrom", signature(object="array"),
	  annotatedDataFrameFromArray)
