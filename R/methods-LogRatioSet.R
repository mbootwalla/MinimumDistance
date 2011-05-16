setMethod("initialize", "LogRatioSet",
	  function(.Object,
		   logRRatio=new("matrix"),...){
		  callNextMethod(.Object,
				 logRRatio=logRRatio,...)

	  })

setAs("BeadStudioSet", "LogRatioSet",
      function(from, to){
	      new("LogRatioSet",
		  logRRatio=logR(from),
		  phenoData=phenoData(from),
		  featureData=featureData(from),
		  protocolData=protocolData(from),
		  annotation=annotation(from))
      })

setMethod("logR", "LogRatioSet", function(object) assayData(object)[["logRRatio"]])
setReplaceMethod("logR", c("LogRatioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "logRRatio", value)
	 })

##setReplaceMethod("logR", c("LogRatioSet", "ffdf"),
##		 function(object, value) {
##			 assayDataElementReplace(object, "logRRatio", value)
##	 })

setMethod("baf", "LogRatioSet",
	  function(object) {
		  assayData(object)[["BAF"]]
	 })

##setReplaceMethod("baf", c("LogRatioSet", "ffdf"),
##		 function(object, value) {
##			 assayDataElementReplace(object, "BAF", value)
##	 })
##setReplaceMethod("baf", c("LogRatioSet", "matrix"),
##		 function(object, value) {
##			 assayDataElementReplace(object, "BAF", value)
##	 })
setReplaceMethod("baf", c("LogRatioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })
