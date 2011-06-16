setGeneric("mindist", function(object) standardGeneric("mindist"))
setGeneric("family", function(object) standardGeneric("family"))
setGeneric("plot", useAsDefault=function(x,y, ...) graphics::plot(x,y,...))
setGeneric("xyplot", useAsDefault=function(x, data, ...) lattice::xyplot(x, data,...))
setGeneric("loglik", function(object) standardGeneric("loglik"))
setGeneric("loglik<-", function(object, value) standardGeneric("loglik<-"))
setGeneric("range.index", function(object) standardGeneric("range.index"))
setGeneric("nMarkers", function(object) standardGeneric("nMarkers"))
setGeneric("coverage", function(object) standardGeneric("coverage"))
setGeneric("state", function(object) standardGeneric("state"))
setGeneric("mindist<-", function(object,value) standardGeneric("mindist<-"))
setGeneric("calculateMindist", function(object, ...) standardGeneric("calculateMindist"))
setGeneric("logR", function(object) standardGeneric("logR"))
setGeneric("logR<-", function(object, value) standardGeneric("logR<-"))
setGeneric("baf<-", function(object, value) standardGeneric("baf<-"))
setGeneric("baf", function(object) standardGeneric("baf"))
setGeneric("mad")
setGeneric("ncol")
setGeneric("mad<-", function(object, value) standardGeneric("mad<-"))
setGeneric("xsegment", function(object, id, ..., verbose=FALSE) standardGeneric("xsegment"))
setGeneric("phenoData2", function(object) standardGeneric("phenoData2"))
setGeneric("phenoData2<-", function(object,value) standardGeneric("phenoData2<-"))
setGeneric("varLabels2", function(object) standardGeneric("varLabels2"))
setGeneric("prune", function(object, ranges, id,
			      lambda=0.05,
			      min.change=0.1,
			      min.coverage=3,
			      scale.exp=0.02,
			      verbose, ...)
	   standardGeneric("prune"))
setGeneric("computeBayesFactor", function(object,
					  ranges,
					  id,
					  states=0:4,
					  baf.sds=c(0.02, 0.03, 0.02),
					  mu.logr=c(-2, -0.5, 0, 0.3, 0.75),
					  log.pi,
					  tau,
					  normal.index=61,
					  a=0.0009,
					  prGtCorrect=0.999,
					  df0=10,
					  verbose=TRUE, ...)
	   standardGeneric("computeBayesFactor"))
setGeneric("todf", function(object, range, frame, ...) standardGeneric("todf"))
setGeneric("offspringNames", function(object) standardGeneric("offspringNames"))
setGeneric("fatherNames", function(object) standardGeneric("fatherNames"))
setGeneric("motherNames", function(object) standardGeneric("motherNames"))
setGeneric("rbind", function(..., deparse.level=1) standardGeneric("rbind"),
           signature = "...")
setGeneric("width", function(x) standardGeneric("width"))

