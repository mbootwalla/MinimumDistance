setMethod("state", signature(object="RangedDataCNV"), function(object) object$state)
setMethod("nMarkers", signature(object="RangedDataCNV"), function(object) object$num.mark)
setMethod("sampleNames", signature(object="RangedDataCNV"), function(object) object$id)

