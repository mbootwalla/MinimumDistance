setOldClass("ffdf")
setOldClass("ff_matrix")
setOldClass("ff_array")

setClass("RangedDataCopyNumber", contains="RangedData",
	 representation("VIRTUAL"))
setClass("RangedDataCNV", contains="RangedDataCopyNumber")
setClass("RangedDataCNVTrios", contains="RangedDataCNV")
## It would be much cleaner to add slots to the RangedData class for
## chrom, id, and num.mark?  This way we can force these slots to be a
## specific class like numeric, integer, etc.
setValidity("RangedDataCNV", function(object){
	all(c("chrom", "id", "num.mark") %in% colnames(object))
})
setClass("RangedDataCBS", contains="RangedDataCNV")
setValidity("RangedDataCBS", function(object) all(c("seg.mean", "start.index", "end.index") %in% colnames(object)))
setClass("RangedDataHMM", contains="RangedDataCNV")
setValidity("RangedDataHMM", function(object) "state" %in% colnames(object))
##setClassUnion("dataFrame", "data.frame")
setClass("DataFrameCNV", contains="data.frame")
setClass("RangedDataCNVList", contains="list")
setMethod("stack", signature(x="RangedDataCNVList"),
	  function(x){
		  x <- lapply(x, function(x){
			  browser()
			  as(x, "RangedData")
		  })
		  rdl <- RangedDataList(x)
		  rd <- stack(rdl)
		  return(rd)
	  })
setClass("RangedDataCBS2", contains="RangedDataCBS")
setValidity("RangedDataCBS2", function(object){
	"state" %in% colnames(object)
})

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
## could include file.ext, cdfname
setClass("SampleSheet", contains="data.frame")
setValidity("SampleSheet", function(object) "Sample.Name" %in% colnames(object))
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

