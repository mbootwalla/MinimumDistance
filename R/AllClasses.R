##setClass("TrioSet", contains="SnpCallSetPlus", representation(isBiparental="logical"))
setClass("LikSet",
	 contains="LogRatioSet",
	 representation(loglik="array",
			range.index="integer"),
	 prototype = prototype(
	 new("VersionedBiobase",
	     versions=c(classVersion("LogRatioSet"), LikSet="1.0.0"))))


setClass("RangedDataCNV", contains="RangedData")
setValidity("RangedDataCNV", function(object){
	all(c("chrom", "id", "num.mark", "state") %in% colnames(object))
})
setClass("RangedDataCNVPlus", contains="RangedDataCNV")
setValidity("RangedDataCNV", function(object){
	all(c("seg.mean", "start.index", "end.index") %in% colnames(object))
})

setClass("MinDistanceSet", contains="MultiSet")

setClassUnion("matrixOrNULL", c("matrix", "NULL"))

setClass("TrioSet", contains="LogRatioSet",
	 representation(phenoData2="array"))

setClass("TrioSet", contains="LogRatioSet",
	 representation(phenoData2="array",
			mindist="matrixOrNULL"),
	 prototype = prototype(
	                       new("VersionedBiobase",
				   versions=c(classVersion("eSet"), TrioSet="0.0.2"))))

setMethod("updateObject", signature(object="TrioSet"),
          function(object, ..., verbose=FALSE) {
		  obj <- tryCatch(callNextMethod(), error=function(e) NULL)
		  if(is.null(obj)){
			  stop("Update was unsucessfull.")
		  }
		  return(obj)
	  })
##  AssayDataElements are T x M arrays
##          T= number features
##          M= 3 (Father, Mother, Offspring)
##
## - the default phenoData can contain information about the trio.
## - phenoData2 is a container for the sample-level metaData
##     -> this must be an T x M x P array,
##        T: number of trios
##        M: father, mother, offspring (3)
##        P: covariates on each sample

## need annotatedDataFrameFromAssayData to call annotatedDataFrameFrom methods for arrays






##options(error=stop)
##new("TrioSet", logRRatio=logRRatio)






##setGeneric("isBiparental", function(object) standardGeneric("isBiparental"))
####setMethod("isBiparental", "TrioSet", function(object) object@isBiparental)
##setMethod("initialize", "TrioSet", function(.Object,
##					    isBiparental, ...){
##	.Object <- callNextMethod(.Object, ...)
##	.Object@isBiparental <- isBiparental
##})


