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
		  if(length(rdList) > 1){
			  rdl <- RangedDataList(rdList)
			  rd <- stack(rdl)
			  rd <- rd[, -ncol(rd)]
		  } else rd <- rdList[[1]]
		  rd <- as(rd, "RangedDataCNV")
		  return(rd)
	  })

setMethod("offspringNames", signature(object="TrioSetList"), function(object) offspringNames(object[[1]]))
setMethod("fatherNames", signature(object="TrioSetList"), function(object) fatherNames(object[[1]]))
setMethod("motherNames", signature(object="TrioSetList"), function(object) motherNames(object[[1]]))

setMethod("offspringNames", signature(object="TrioSet"), function(object){
	phenoData2(object)[, "CIDR_Name", "O"]
})
setMethod("fatherNames", signature(object="TrioSet"), function(object){
	phenoData2(object)[, "CIDR_Name", "F"]
})

setMethod("motherNames", signature(object="TrioSet"), function(object){
	phenoData2(object)[, "CIDR_Name", "M"]
})


