##setMethod("baf", "MultiSet", function(object) assayData(object)[["BAF"]])
##setMethod("logR", "MultiSet", function(object) assayData(object)[["logR"]])
##
##setMethod("baf", "LogRatioSet", function(object) assayData(object)[["BAF"]])
##setMethod("logR", "LogRatioSet", function(object) assayData(object)[["logRRatio"]])
##
##
##setMethod("logR.M", "MultiSet", function(object) assayData(object)[["logR.M"]])
##setMethod("logR.F", "MultiSet", function(object) assayData(object)[["logR.F"]])
##setMethod("logR.O", "MultiSet", function(object) assayData(object)[["logR.O"]])
##setMethod("baf.M", "MultiSet", function(object) assayData(object)[["baf.M"]])
##setMethod("baf.F", "MultiSet", function(object) assayData(object)[["baf.F"]])
##setMethod("baf.O", "MultiSet", function(object) assayData(object)[["baf.O"]])
##setMethod("mindist", "MultiSet", function(object) assayData(object)[["mindist"]])
##
##setMethod("copyNumber", "ExpressionSet", function(object) exprs(object))
##setMethod("position", "MultiSet", function(object) fData(object)$position)
##setMethod("chromosome", "MultiSet", function(object) fData(object)$chromosome)
##
##
##setAs("MultiSet", "LogRatioSet",
##      function(from, to){
##	      phenodata=phenoData(from)
##	      sampleNames(phenodata)=phenodata$Sample.Name
##	      prD <- protocolData(from)
##	      sampleNames(prD) <- sampleNames(phenodata)
##	      lR <- assayData(from)[["logR"]]
##	      ix <- match(colnames(lR), sampleNames(phenodata))
##	      phenodata <- phenodata[ix, ]
##	      prD <- prD[ix, ]
##	      new("LogRatioSet",
##		  logRRatio=assayData(from)[["logR"]],
##		  phenoData=phenodata,
##		  protocolData=prD,
##		  featureData=featureData(bsSet),
##		  annotation=annotation(bsSet))
##      })
##
##setAs("MultiSet", "data.frame",
##	  function(from){
##		  data.frame(logR=c(logR.F(from),
##			     logR.M(from),
##			     logR.O(from),
##			     -mindist(from)),
##			     baf=c(baf.F(from),
##			     baf.M(from),
##			     baf.O(from),
##			     rep(NA, nrow(from))),
##			     x=rep(position(from)/1e6, 4),
##			     subject=factor(rep(1:4, each=nrow(from)),
##			                    labels=c("Father", "Mother", "Offspring", "distance"),
##					    ordered=TRUE))
##	  })


setMethod("family", "eSet", function(object) substr(sampleNames(object), 1, 5))
