setMethod("baf", "MultiSet", function(object) assayData(object)[["BAF"]])
setMethod("logR", "MultiSet", function(object) assayData(object)[["logR"]])
