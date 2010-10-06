setMethod("baf", "MultiSet", function(object) assayData(object)[["BAF"]])
setMethod("logR", "MultiSet", function(object) assayData(object)[["logR"]])

setMethod("baf", "LogRatioSet", function(object) assayData(object)[["BAF"]])
setMethod("logR", "LogRatioSet", function(object) assayData(object)[["logRRatio"]])
