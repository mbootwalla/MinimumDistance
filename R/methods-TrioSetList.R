setMethod("xsegment", signature(object="TrioSetList"),
	  function(object, id, ..., verbose=FALSE){
		  dfl <- lapply(object, xsegment, id=id, ..., verbose=verbose)
		  df <- do.call("rbind", dfl)
		  return(df)
	  })
setMethod("sampleNames", signature(object="TrioSetList"),
	  function(object) sampleNames(object[[1]]))
setMethod("ncol", signature(x="TrioSetList"),
	  function(x) ncol(x[[1]]))
setMethod("prune", signature(object="TrioSetList", ranges="RangedDataCNV"),
	  function(object, ranges, id, lambda, min.change, min.coverage,
		   scale.exp, verbose, ...){
		  rdList <- lapply(object, prune, ranges=ranges,
				   id=id,
				   lambda=lambda,
				   min.change=min.change,
				   min.coverage=min.coverage,
				   scale.exp=scale.exp,
				   verbose=verbose, ...)
	  })
