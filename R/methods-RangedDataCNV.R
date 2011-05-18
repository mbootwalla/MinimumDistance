setMethod("state", signature(object="RangedDataCNV"), function(object) object$state)
setMethod("nMarkers", signature(object="RangedDataCNV"), function(object) object$num.mark)
setMethod("sampleNames", signature(object="RangedDataCNV"), function(object) object$id)
setMethod("RangedDataCNV", signature(ranges="IRanges"),
	  function(ranges=IRanges(), ...,
		   space=NULL,
		   universe=NULL){
		  nms <- names(list(...))
		  stopifnot(c("chrom", "id", "num.mark", "state", "seg.mean", "start.index", "end.index") %in% nms)
		  rd <- RangedData(ranges=ranges, ..., space=space, universe=universe)
		  rd2 <- as(rd, "RangedDataCNV")
		  return(rd2)
	  })


