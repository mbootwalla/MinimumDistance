setMethod("baf", "MultiSet", function(object) assayData(object)[["BAF"]])
setMethod("logR", "MultiSet", function(object) assayData(object)[["logR"]])

setMethod("baf", "LogRatioSet", function(object) assayData(object)[["BAF"]])
setMethod("logR", "LogRatioSet", function(object) assayData(object)[["logRRatio"]])


setMethod("logR.M", "MultiSet", function(object) assayData(object)[["logR.M"]])
setMethod("logR.F", "MultiSet", function(object) assayData(object)[["logR.F"]])
setMethod("logR.O", "MultiSet", function(object) assayData(object)[["logR.O"]])
setMethod("baf.M", "MultiSet", function(object) assayData(object)[["baf.M"]])
setMethod("baf.F", "MultiSet", function(object) assayData(object)[["baf.F"]])
setMethod("baf.O", "MultiSet", function(object) assayData(object)[["baf.O"]])
setMethod("mindist", "MultiSet", function(object) assayData(object)[["mindist"]])

setMethod("copyNumber", "ExpressionSet", function(object) exprs(object))
setMethod("position", "MultiSet", function(object) fData(object)$position)
setMethod("chromosome", "MultiSet", function(object) fData(object)$chromosome)
