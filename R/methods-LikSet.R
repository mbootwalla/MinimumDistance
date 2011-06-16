setMethod("range.index", signature(object="LikSet"), function(object) fData(object)$range.index)
assayDataStorageMode <- Biobase:::assayDataStorageMode
assayDataEnvLock <- Biobase:::assayDataEnvLock
setMethod("[", "LikSet", function(x, i, j, ..., drop=FALSE){
	## why does callNextMethod seem to call itself and not the eset method
	##x <- callNextMethod(x, i, j, ..., drop=drop)
	if (missing(drop))
		drop <- FALSE
	if (missing(i) && missing(j)) {
		if (length(list(...))!=0)
			stop("specify genes or samples to subset; use '",
			     substitute(x), "$", names(list(...))[[1]],
			     "' to access phenoData variables")
		return(x)
	}
	if (!missing(j)) {
		phenoData(x) <- phenoData(x)[j,, ..., drop = drop]
		protocolData(x) <- protocolData(x)[j,, ..., drop = drop]
	}
	if (!missing(i))
		featureData(x) <- featureData(x)[i,,..., drop=drop]
	## assayData; implemented here to avoid function call
	orig <- assayData(x)
	storage.mode <- assayDataStorageMode(orig)
	assayData(x) <-
		switch(storage.mode,
		       environment =,
		       lockedEnvironment = {
			       aData <- new.env(parent=emptyenv())
			       if (missing(i))                     # j must be present
				       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][, j, ..., drop = drop]
			       else {                              # j may or may not be present
				       if (missing(j))
					       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i,, ..., drop = drop]
				       else
					       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i, j, ..., drop = drop]
			       }
			       if ("lockedEnvironment" == storage.mode) assayDataEnvLock(aData)
			       aData
		       },
		       list = {
			       if (missing(i))                     # j must be present
				       lapply(orig, function(obj) obj[, j, ..., drop = drop])
			       else {                              # j may or may not be present
				       if (missing(j))
					       lapply(orig, function(obj) obj[i,, ..., drop = drop])
				       else
					       lapply(orig, function(obj) obj[i, j, ..., drop = drop])
			       }
		       })
	loglik(x) <- loglik(x)[ , i, j, , ..., drop=drop]
	x
})

setMethod("loglik", signature(object="LikSet"), function(object) object@loglik)
setReplaceMethod("loglik", signature(object="LikSet", value="array"), function(object,value){
	object@loglik <- value
	object
})
setReplaceMethod("loglik", signature(object="LikSet", value="numeric"), function(object,value){
	object@loglik <- value
	object
})
