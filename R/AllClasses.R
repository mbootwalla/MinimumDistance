setOldClass("ffdf")
setOldClass("ff_matrix")
setOldClass("ff_array")
setClass("RangedDataCopyNumber", contains="RangedData",
	 representation("VIRTUAL"))
setClass("RangedDataCNV", contains="RangedDataCopyNumber")
setValidity("RangedDataCNV", function(object){
	##all(c("chrom", "id", "num.mark",  "start.index", "end.index") %in% colnames(object))
	is(object, "RangedData")
})
setClass("RangedDataCBS", contains="RangedDataCNV")
setValidity("RangedDataCBS", function(object){
	"seg.mean" %in% colnames(object)
})
setClass("RangedDataHMM", contains="RangedDataCNV")
setValidity("RangedDataHMM", function(object){
	"state" %in% colnames(object)
})
##setMethod("initialize", "RangedDataCNV",
##	  function(.Object,
##		   ranges, values, ...){
##		  .Object <- callNextMethod(.Object, ranges=ranges, values=values, ...)
##	  })
RangedDataCNV <- function(ranges=IRanges(),
			  start,
			  end,
			  chromosome,
			  coverage,
			  sampleId,
			  startIndexInChromosome,
			  endIndexInChromosome,
			  ...){
	if(!missing(end) && !missing(start))
		ranges <- IRanges(start, end)
	if(missing(chromosome))
		chromosome <- vector("integer", length(ranges))
	if(missing(coverage))
		coverage <- vector("integer", length(ranges))
	if(missing(sampleId))
		sampleId <- vector("character", length(ranges))
	if(missing(startIndexInChromosome))
		startIndexInChromosome <- vector("integer", length(ranges))
	if(missing(endIndexInChromosome))
		endIndexInChromosome <- vector("integer", length(ranges))
	rd <- RangedData(ranges,
			 chrom=chromosome,
			 num.mark=coverage,
			 id=sampleId,
			 start.index=startIndexInChromosome,
			 end.index=endIndexInChromosome, ...)##, ...)
	new("RangedDataCNV", ranges=ranges(rd), values=values(rd))
}

RangedDataCBS <- function(ranges=IRanges(),
			  seg.mean=vector("numeric", length(ranges)), ...){
	rd <- RangedDataCNV(ranges=ranges, seg.mean=seg.mean, ...)
	new("RangedDataCBS", ranges=ranges(rd), values=values(rd))
}
RangedDataHMM <- function(ranges=IRanges(),
			  state=vector("integer", length(ranges)), ...){
	rd <- RangedDataCNV(ranges=ranges, state=state, ...)
	new("RangedDataHMM", ranges=ranges(rd), values=values(rd))
}

setClass("MinDistanceSet", contains="MultiSet")

setClassUnion("matrixOrNULL", c("matrix", "NULL", "ff_matrix"))
setClassUnion("arrayOrNULL", c("array", "NULL"))

setClass("LogRatioSet", contains="eSet")
setClass("BeadStudioSet", contains="eSet")

setClass("LikSet",
	 contains="LogRatioSet",
	 representation(loglik="array",
			range.index="integer"),
	 prototype = prototype(
	 new("VersionedBiobase",
	     versions=c(classVersion("LogRatioSet"), LikSet="1.0.0"))))

##setClass("TrioSet", contains="LogRatioSet",
##	 representation(phenoData2="array"))
##
##setClass("TrioSet", contains="LogRatioSet",
##	 representation(phenoData2="array",
##			mindist="matrixOrNULL"),
##	 prototype = prototype(
##	                       new("VersionedBiobase",
##				   versions=c(classVersion("eSet"), TrioSet="0.0.2"))))
##
##
##setClass("TrioSet", contains="LogRatioSet",
##	 representation(phenoData2="arrayOrNULL",
##			mindist="matrixOrNULL",
##			mad="matrix"),
##	 prototype = prototype(
##	                       new("VersionedBiobase",
##				   versions=c(classVersion("LogRatioSet"), TrioSet="0.0.3"))))

setClass("TrioSet", contains="LogRatioSet",
	 representation(phenoData2="arrayOrNULL",
			mindist="matrixOrNULL",
			mad="matrix"),
	 prototype = prototype(
	                       new("VersionedBiobase",
				   versions=c(classVersion("LogRatioSet"), TrioSet="0.0.4"))))

## should we add a slot for trioNames
##  -- would be a R x 3 matrix, where R is the number of trios
##  -- R must be equal to the number of columns of the assayData arrays
##  -- '[' method for the trioNames slot
##  -- 'show' method for trioNames slot
##  -- add trioNames accessor / replacement method
##  -- modify offspringNames, fmoNames, ... to access trioNames slot
## Other potential slots:
##  -- dna source
##  -- batch / plate
setClass("TrioSetList", contains="list")
##setClass("TrioSetList",
##	 representation(trioSets="list",
##			phenoData="AnnotatedDataFrame",
##			mad="matrix",
##			dnaSource="matrixOrNull",
##			trioNames="matrixOrNull",
##			batch="matrixOrNull"))


##setClass("TrioSet", contains="BeadStudioSet",
##	 representation(phenoData2="array",
##			mindist="matrixOrNULL"),
##	 prototype = prototype(
##	                       new("VersionedBiobase",
##				   versions=c(classVersion("eSet"), TrioSet="0.0.3"))))

setMethod("updateObject", signature(object="TrioSet"),
          function(object, ..., verbose=FALSE) {
		  obj <- tryCatch(callNextMethod(), error=function(e) NULL)
		  if(is.null(obj)){
			  stop("updateObject failed")
##			  md <- tryCatch(mindist(object), error=function(e) NULL)
##			  if(is.null(md)){
##				  object <- new("TrioSet",
##					     logRRatio=logR(object),
##					     BAF=baf(object),
##					     phenoData=phenoData(object),
##					     phenoData2=object@phenoData2,
##					     experimentData=experimentData(object),
##					     featureData=featureData(object),
##					     protocolData=protocolData(object),
##					     mindist=NULL,
##					     annotation=annotation(object))
##				  return(object)
##			  } else {
##
##			  }
##			  mads <- tryCatch(mad(object), error=function(e) NULL)
##			  if(is.null(mads)){
##				  callNextMethod(mad=array())
##			  } else
##				  object <- new("TrioSet",
##					     logRRatio=logR(object),
##					     BAF=baf(object),
##					     phenoData=phenoData(object),
##					     phenoData2=object@phenoData2,
##					     experimentData=experimentData(object),
##					     featureData=featureData(object),
##					     protocolData=protocolData(object),
##					     mindist=mindist(object),
##					     annotation=annotation(object))
##			  }
		  }
		  return(object)
	  })

##setClassUnion("dataFrame", "data.frame")
setClass("DataFrameCNV", contains="data.frame")
##setMethod("initialize", signature(.Object="DataFrameCNV"),
##	  function(.Object, ...){
##		  .Object <- callNextMethod(.Object, ...)
##	  })
##DataFrameCNV <- function(...) data.frame()



##setClass("DataFrameCNV", representation(row.names="character",
##					names="character"),
##	 contains="list")


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


