setMethod("initialize", signature(.Object="TrioSet"),
	  function(.Object, phenoData2, ...){
		  .Object <- callNextMethod(.Object, ...)
		  phenoData2 <- array(NA, dim=c(ncol(.Object), 3, ncol(phenoData(.Object))),
				      dimnames=list(sampleNames(.Object),
				      c("F", "M", "O"), NULL))
		  .Object@phenoData2 <- phenoData2
		  .Object
})

setMethod("show", signature(object="TrioSet"), function(object){
	cat(class( object ), " (storageMode: ", storageMode(object), ")\n", sep="")
	adim <- dim(object)
	if (length(adim)>1)
		cat("assayData:",
		    if (length(adim)>1)
		    paste(adim[[1]], "Features,",
			  adim[[2]], "Trios, ",
			  adim[[3]], "Family members (Father, Mother, Offspring)") else NULL,
		    "\n")
	cat("  element names:",
	    paste(assayDataElementNames(object), collapse=", "), "\n")
	Biobase:::.showAnnotatedDataFrame(protocolData(object),
				labels=list(object="protocolData"))
	Biobase:::.showAnnotatedDataFrame(phenoData(object),
				labels=list(object="phenoData"))
	Biobase:::.showAnnotatedDataFrame(featureData(object),
				labels=list(
				object="featureData",
				sampleNames="featureNames",
				varLabels="fvarLabels",
				varMetadata="fvarMetadata"))
	cat("experimentData: use 'experimentData(object)'\n")
	pmids <- pubMedIds(object)
	if (length(pmids) > 0 && all(pmids != ""))
		cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
	cat("Annotation:", annotation(object), "\n")
	cat("mindist:", class(mindist(object)), "\n")
	cat("phenoData2:", dimnames(object@phenoData2)[[3]], "\n")
})

setReplaceMethod("sampleNames", signature(object="TrioSet"),
		 function(object,value){
			 object <- callNextMethod(object, value)
			 row.names(object@phenoData2) <- value
			 object
		 })

setMethod("mindist", "TrioSet", function(object) object@mindist)
setReplaceMethod("mindist", signature(object="TrioSet", value="ANY"),
		 function(object, value){
			 object@mindist <- value
			 return(object)
		 })

##setReplaceMethod("trioNames", signature(object="TrioSet"),
##		 function(object,value){
##			 object <- callNextMethod(object, value)
##			 row.names(object@phenoData2) <- value
##			 object
##		 })

setMethod("dim", "TrioSet", function(x) assayDataDim(assayData(x)))


## Fix strange behavor for [ with ff_arrays.
setMethod("[", signature(x="ff_array"), function(x, i, j, ..., drop=FALSE){
	if(missing(drop)) drop <- TRUE
	if(length(list(...)) > 0){
		k <- list(...)[[1]]
		if(is(k, "character")) {
			k <- match(dimnames(x)[[3]])
			stopifnot(length(k)>0)
		}
	} else k <- NULL
	if(is.null(k)){
		if(!missing(i) & missing(j)){
			x <- x[i, 1:ncol(x), 1:3, ..., drop=drop]
		}
		if(!missing(i) & !missing(j)){
			x <- x[1:nrow(x), j, 1:3, drop=drop]
		}
		if(missing(i) & !missing(j)){
			x <- x[1:nrow(x), j, 1:3, drop=drop]
		}
	} else {
		if(!missing(i) & missing(j)){
			x <- x[i, 1:ncol(x), k, drop=drop]
		}
		if(!missing(i) & !missing(j)){
			x <- x[i, j, k, drop=drop]
		}
		if(missing(i) & !missing(j)){
			x <- x[1:nrow(x), j, k, drop=drop]
		}
	}
	return(x)
})

setMethod("[", "TrioSet", function(x, i, j, ..., drop = FALSE) {
	if (missing(drop))
		  drop <- FALSE
	if(length(list(...)) > 0){
		k <- list(...)[[1]]
	} else k <- NULL
	if (missing(i) && missing(j) && is.null(k)) {
		if (length(list(...))!=0)
			stop("specify features, trios, or samples to subset; use '",
			     substitute(x), "$", names(list(...))[[1]],
			     "' to access phenoData variables")
		return(x)
	}
	if (!missing(j)) {
		phenoData(x) <- phenoData(x)[j,, ..., drop = drop]
		protocolData(x) <- protocolData(x)[j,, ..., drop = drop]
	}
	if (!missing(i))
		featureData(x) <- featureData(x)[i, ,..., drop=drop]
	if (!is.null(k) && !missing(j)){
		x@phenoData2 <- x@phenoData2[j, k, , drop=drop]
	}
	if(!missing(j) && is.null(k)){
		x@phenoData2 <- x@phenoData2[j, , , drop=drop]
	}
	if(missing(j) && !is.null(k)){
		x@phenoData2 <- x@phenoData2[, k, , drop=drop]
	}
	if(!missing(i) & missing(j)){
		mindist(x) <- mindist(x)[i, ...,drop=drop]
	}
	if(!missing(i) & !missing(j)){
		mindist(x) <- mindist(x)[i, j, drop=drop]
	}
	if(missing(i) & !missing(j)){
		mindist(x) <- mindist(x)[, j, drop=drop]
	}
	## assayData; implemented here to avoid function call
	orig <- assayData(x)
	storage.mode <- Biobase:::assayDataStorageMode(orig)
	if(is.null(k)){
		assayData(x) <-
			switch(storage.mode,
			       environment =,
			       lockedEnvironment = {
				       aData <- new.env(parent=emptyenv())
				       if (missing(i))                     # j must be present
					       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][, j, , ..., drop = drop]
				       else {                              # j may or may not be present
					       if (missing(j))
						       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i,, , ..., drop = drop]
					       else
						       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i, j, , ..., drop = drop]
				       }
				       if ("lockedEnvironment" == storage.mode) Biobase:::assayDataEnvLock(aData)
				       aData
			       },
			       list = {
				       if (missing(i))                     # j must be present
					       lapply(orig, function(obj) obj[, j, , ..., drop = drop])
				       else {                              # j may or may not be present
					       if (missing(j))
						       lapply(orig, function(obj) obj[i,, , ..., drop = drop])
					       else
						       lapply(orig, function(obj) obj[i, j, , ..., drop = drop])
				       }
			       })
		return(x)
	} else{ ## not missing k
		assayData(x) <-
			switch(storage.mode,
			       environment =,
			       lockedEnvironment = {
				       aData <- new.env(parent=emptyenv())
				       if (missing(i))                     # j must be present
					       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][, j, k, drop = drop]
				       else {                              # j may or may not be present
					       if (missing(j))
						       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i,, k, drop = drop]
					       else
						       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i, j, k, drop = drop]
				       }
				       if ("lockedEnvironment" == storage.mode) Biobase:::assayDataEnvLock(aData)
				       aData
			       },
			       list = {
				       if (missing(i))                     # j must be present
					       lapply(orig, function(obj) obj[, j, k, ..., drop = drop])
				       else {                              # j may or may not be present
					       if (missing(j))
						       lapply(orig, function(obj) obj[i,, k, ..., drop = drop])
					       else
						       lapply(orig, function(obj) obj[i, j, k, ..., drop = drop])
				       }
			       })
		return(x)
	}
})



setReplaceMethod("logR", signature(object="TrioSet", value="ANY"),
		 function(object, value){
			 assayDataElementReplace(object, "logRRatio", value)
		 })
setReplaceMethod("baf", signature(object="TrioSet", value="ANY"),
		 function(object, value){
			 assayDataElementReplace(object, "BAF", value)
		 })

fullId <- function(object) object@phenoData2[, "Sample.Name", ]

setMethod("calculateMindist", signature(object="TrioSet"), function(object, ...){
	sns <- sampleNames(object)
	##sns <- ssampleNames(bsSet)
	##sns.all <- sampleNames(bsSet)
	##trios <- completeTrios(bsSet)
	##father <- match(trios[, "F"], sns)
	##mother <- match(trios[, "M"], sns)
	##offspring <- match(trios[, "O"], sns)
	##stopifnot(identical(sns[offspring], ssampleNames(minDistanceSet)))
	invisible(open(logR(object)))
	invisible(open(mindist(object)))
	object$MAD <- NA
	##invisible(open(copyNumber(minDistanceSet)))
	##min.resid <- rep(NA, nrow(bsSet))
	for(j in seq(length=ncol(object))){
		##for(j in seq_along(father)){ ## column by column so as not to swamp RAM
		if(j %% 100 == 0) cat(".")
##		f <- father[j]
##		m <- mother[j]
##		o <- offspring[j]
##		logR.f <- logR(bsSet)[, f]
##		logR.m <- logR(bsSet)[, m]
##		logR.o <- logR(bsSet)[, o]
		lr <- logR(object)[, j, ]
		d1 <- lr[, "F"] - lr[, "O"]
		d2 <- lr[, "M"] - lr[, "O"]
		I <- as.numeric(abs(d1) <= abs(d2))
		md <- I*d1 + (1-I)*d2
		mindist(object)[, j] <- md
		##copyNumber(minDistanceSet)[,j] <- md
		object$MAD[j] <- mad(md, na.rm=TRUE)
		##minDistanceSet$MAD[j] <- mad(md, na.rm=TRUE)
	}
	close(mindist(object))
	close(logR(object))
	##mindist[index, j] <- d2[-index]
##	for(j in seq_along(sns)){
##		offspring.name <- sns[j]
##		family <- substr(offspring.name, 1, 5)
##		father.name <- paste(family, "03", sep="_")
##		father.index <- match(father.name, sns.all)
##		mother.name <- paste(family, "02", sep="_")
##		mother.index <- match(mother.name, sns.all)
##		offspring.index <- match(offspring.name, sns.all)
##		lR <- as.matrix(logR(bsSet)[, c(mother.index, father.index)])
##		resid <- lR - logR(bsSet)[, offspring.index]
##		na.index <- which(rowSums(is.na(resid)) > 0)
##		##resid[is.na(resid)] <- 0
##		sign.resid <- sign(resid)
##		min.resid[-na.index] <- rowMin(abs(resid[-na.index, ]))
##		## now give the appropriate sign
##		column.index <- ifelse(abs(resid[,1]) < abs(resid[, 2]), 1, 2)
##		col2 <- which(column.index == 2)
##		col1 <- which(column.index == 1)
##		min.resid[col2] <- min.resid[col2] * sign.resid[col2, 2]
##		min.resid[col1] <- min.resid[col1] * sign.resid[col1, 1]
##		##min.resid <- rowMin(resid)
##		copyNumber(minDistanceSet)[, j] <- min.resid
##		if(j %% 100 == 0) cat(j, " ")
##		minDistanceSet$MAD[j] <- mad(min.resid, na.rm=TRUE)
##	}
	return(object)
})

setMethod("mad", signature(x="TrioSet"), function(x) x@mad)

setReplaceMethod("mad", signature(object="TrioSet", value="ANY"),
	  function(object, value){
		  object@mad <- value
		  return(object)
	  })
