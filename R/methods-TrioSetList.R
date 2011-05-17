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
