setGeneric("addSampleSheet", function(object) standardGeneric("addSampleSheet"))
setGeneric("axis.limits", function(object) standardGeneric("axis.limits"))
setGeneric("biparentalHmm", function(object, ...) standardGeneric("biparentalHmm"))
setGeneric("computeBpiEmission", function(object, ...) standardGeneric("computeBpiEmission"))
setGeneric("gtConfidence", function(object) standardGeneric("gtConfidence"))
setGeneric("inheritance", function(object) standardGeneric("inheritance"))
setGeneric("informative", function(object) standardGeneric("informative"))
setGeneric("mendelian", function(object) standardGeneric("mendelian"))
setGeneric("isBiparental", function(object, ...) standardGeneric("isBiparental"))
##setGeneric("callDeletion", function(object, ...) standardGeneric(object))

##setGeneric("[", function(x, i, j, k, ..., drop=FALSE) standardGeneric("["))

setGeneric("logR.M", function(object) standardGeneric("logR.M"))
setGeneric("logR.F", function(object) standardGeneric("logR.F"))
setGeneric("logR.O", function(object) standardGeneric("logR.O"))
setGeneric("baf.M", function(object) standardGeneric("baf.M"))
setGeneric("baf.F", function(object) standardGeneric("baf.F"))
setGeneric("baf.O", function(object) standardGeneric("baf.O"))
setGeneric("mindist", function(object) standardGeneric("mindist"))
setGeneric("family", function(object) standardGeneric("family"))

setGeneric("plot", useAsDefault=function(x,y, ...) graphics::plot(x,y,...))
setGeneric("xyplot", useAsDefault=function(x, data, ...) lattice::xyplot(x, data,...))
##setGeneric("plot", function(object, ...) standardGeneric("plot"))
##setGeneric("plot")

setGeneric("loglik", function(object) standardGeneric("loglik"))
setGeneric("loglik<-", function(object, value) standardGeneric("loglik<-"))
setGeneric("range.index", function(object) standardGeneric("range.index"))

setGeneric("nMarkers", function(object) standardGeneric("nMarkers"))
setGeneric("state", function(object) standardGeneric("state"))
setGeneric("mindist<-", function(object,value) standardGeneric("mindist<-"))
setGeneric("calculateMindist", function(object, ...) standardGeneric("calculateMindist"))

##setGeneric("updateObject", function(object, ...) standardGeneric("updateObject"))
setGeneric("logR", function(object) standardGeneric("logR"))
setGeneric("logR<-", function(object, value) standardGeneric("logR<-"))
setGeneric("baf<-", function(object, value) standardGeneric("baf<-"))
setGeneric("baf", function(object) standardGeneric("baf"))

setGeneric("mad")
setGeneric("ncol")
setGeneric("mad<-", function(object, value) standardGeneric("mad<-"))
##setGeneric("segment")
setGeneric("xsegment", function(object, id, ..., verbose=FALSE) standardGeneric("xsegment"))

setGeneric("phenoData2", function(object) standardGeneric("phenoData2"))
setGeneric("phenoData2<-", function(object,value) standardGeneric("phenoData2<-"))
setGeneric("varLabels2", function(object) standardGeneric("varLabels2"))

setGeneric("RangedDataCNV", function(ranges=IRanges(), ..., space=NULL, universe=NULL)
	   standardGeneric("RangedDataCNV"))


setGeneric("prune", function(object, ranges, id,
			      lambda=0.05,
			      min.change=0.1,
			      min.coverage=3,
			      scale.exp=0.02,
			      verbose, ...)
	   standardGeneric("prune"))


setGeneric("todf", function(object, col=1:3, ...) standardGeneric("todf"))
