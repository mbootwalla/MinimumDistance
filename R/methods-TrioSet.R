setMethod("initialize", signature(.Object="TrioSet"),
	  function(.Object, phenoData2, ...){
		  .Object <- callNextMethod()
		  if(!"logRRatio" %in% names(list(...))){
			  .Object <- assayDataElementReplace(.Object, "logRRatio", array(NA, dim=c(0,0,0)))
		  }
		  if("phenoData2" %in% names(list(...))){
			  phenoData2 <- list(...)[["phenoData2"]]
		  } else phenoData2 <- NULL
		  .Object@phenoData2 <- phenoData2
		  .Object@mad <- matrix(NA, ncol(.Object), 3)
		  dimnames(.Object@mad) <- list(sampleNames(.Object), c("F", "M", "O"))
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
	cat("phenoData2:", class(object@phenoData2), "\n")
	cat("mad:", class(mad(object)), "\n")
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
		mad(x) <- mad(x)[j,,...,drop=drop]
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

setMethod("phenoData2", signature(object="TrioSet"), function(object) object@phenoData2)
setReplaceMethod("phenoData2", signature(object="TrioSet", value="ANY"), function(object, value){
	object@phenoData2 <- value
	object
})
setMethod("varLabels2", signature(object="TrioSet"), function(object) colnames(phenoData2(object)))


setMethod("xsegment", signature(object="TrioSet"),
	  function(object, id, ..., verbose=FALSE){
		  ## needs to be ordered
		  ix <- order(chromosome(object), position(object))
		  stopifnot(all(diff(ix) > 0))
		  ##
		  ##
		  ##dfl <- vector("list", 22) ## data.frame list
		  if(missing(id)) id <- sampleNames(object)
		  sample.index <- match(id, sampleNames(object))
		  stopifnot(length(sample.index) > 0)
		  open(mindist(object))
		  fns <- featureNames(object)
		  ##
		  ##
		  ##marker.index.list <- split(seq(length=nrow(minDistanceSet)), chromosome(minDistanceSet))
		  ##for(CHR in 1:22){
		  marker.index <- seq(length=nrow(object))[!duplicated(position(object))]
		  pos <- position(object)[marker.index]
		  chrom <- chromosome(object)[marker.index]
		  CN <- mindist(object)[marker.index, sample.index, drop=FALSE]
		  CN <- matrix(as.numeric(CN), nrow(CN), ncol(CN))
		  dimnames(CN) <- list(featureNames(object)[marker.index], sampleNames(object)[sample.index])
		  arm <- getChromosomeArm(chrom, pos)
		  index.list <- split(seq_along(marker.index), arm)
		  md.segs <- list()
		  if(verbose) message("Running CBS by chromosome arm")
		  for(i in seq_along(index.list)){
			  j <- index.list[[i]]
			  CNA.object <- CNA(genomdat=CN[j, , drop=FALSE],
					    chrom=chrom[j],
					    maploc=pos[j],
					    data.type="logratio",
					    sampleid=id)
			  smu.object <- smooth.CNA(CNA.object)
			  tmp <- segment(smu.object, verbose=0, ...)
			  df <- tmp$output
			  sr <- tmp$segRows
			  ##df <- cbind(tmp$output, tmp$segRows)
			  ##md.segs[[i]] <-
			  firstMarker <- rownames(CNA.object)[sr$startRow]
			  endMarker <- rownames(CNA.object)[sr$endRow]
			  df$start.index <- match(firstMarker, fns)
			  df$end.index <- match(endMarker, fns)
			  md.segs[[i]] <- df
		  }
		  if(length(md.segs) > 1){
			  md.segs <- do.call("rbind", md.segs)
		  } else md.segs=md.segs[[1]]
		  close(mindist(object))
		  return(md.segs)
})

setMethod("computeBayesFactor", signature(object="TrioSet"),
	  function(object, ranges, id,
		   states,
		   baf.sds,
		   mu.logr,
		   log.pi,
		   tau,
		   normal.index,
		   a,
		   prGtCorrect,
		   df0,
		   verbose, ...){
	stopifnot(!missing(tau))
	stopifnot(!missing(log.pi))
	stopifnot(!missing(object))
	CHR <- unique(chromosome(object))
	ranges <- ranges[chromosome(ranges) == CHR, ]
	if(missing(id)) id <- unique(ranges$id) else stopifnot(id %in% unique(ranges$id))
	ranges <- ranges[ranges$id %in% id, ]
	ranges$bayes.factor <- NA
	ranges$argmax <- NA
	ranges$DN <- NA
	for(i in seq_along(id)){
		if(verbose)
			message("\t bayes factors for sample ", i, " of ", length(id))
		this.id <- id[i]
		if(verbose){
			if(i %% 100 == 0)
				message("   sample ", this.id, " (", i, "/", length(id), ")")
		}
		j <- which(ranges$id == this.id)
		rd <- joint4(trioSet=object,
			     ranges=ranges[j, ],
			     states=states,
			     baf.sds=baf.sds,
			     mu.logr=mu.logr,
			     log.pi=log.pi,
			     tau=tau,
			     normal.index=normal.index,
			     a=a,
			     prGtCorrect=prGtCorrect,
			     df0=df0,
			     verbose=verbose)##, F=F, M=M, O=O)
		ranges$bayes.factor[j] <- rd$bayes.factor
		ranges$argmax[j] <- rd$argmax
		ranges$DN[j] <- rd$DN
	}
	ranges
})

setMethod("todf", signature(object="TrioSet", range="RangedData"),
	  function(object, range, frame, ...){
	stopifnot(nrow(range) == 1)
	if(missing(frame)){
		w <- width(range)
		frame <- w/0.05  * 1/2
	}
	marker.index <- Beaty:::featuresInRange(object, range, FRAME=frame)
	id <- range$id
	sample.index <- match(id, offspringNames(object))
	stopifnot(length(sample.index)==1)
	open(baf(object))
	open(logR(object))
	open(mindist(object))
	b <- baf(object)[marker.index, sample.index, ]
	r <- logR(object)[marker.index, sample.index, ]
	md <- mindist(object)[marker.index, sample.index]
	close(baf(object))
	close(logR(object))
	close(mindist(object))
	id <- matrix(c("father", "mother", "offspring"), nrow(b), ncol(b), byrow=TRUE)
	empty <- rep(NA, length(md))
	## A trick to add an extra panel for genes and cnv
	##df <- rbind(df, list(as.integer(NA), as.numeric(NA), as.numeric(NA), as.factor("genes")))
	## The NA's are to create extra panels (when needed for lattice plotting)
##	b <- c(as.numeric(b), empty, 0, 0)
##	r <- c(as.numeric(r), md, 0, 0)
##	x <- c(rep(position(object)[marker.index], 4), 0, 0)
	id <- c(as.character(id), rep("min dist",length(md)))##, c("genes", "CNV"))
	b <- c(as.numeric(b), empty)
	r <- c(as.numeric(r), md)
	x <- rep(position(object)[marker.index], 4)/1e6
	df <- data.frame(x=x, b=b, r=r, id=id)

	df2 <- data.frame(id=c(as.character(df$id), "genes", "CNV"),
			  b=c(df$b, NA, NA),
			  r=c(df$r, NA, NA),
			  x=c(df$x, NA, NA))
	df2$id <- factor(df2$id, levels=c("father", "mother", "offspring", "min dist", "genes", "CNV"), ordered=TRUE)
	return(df2)
})

setMethod("prune", signature(object="TrioSet", ranges="RangedDataCNV"),
	  function(object, ranges, ...){
		  CHR <- unique(chromosome(object))
		  if(verbose) message("Pruning chromosome ", CHR)
		  if(missing(id)) id <- unique(sampleNames(ranges))
		  index <- which(chromosome(ranges) == CHR & sampleNames(ranges) %in% id)
		  ranges <- ranges[index, ]
		  rdList <- vector("list", length(unique(id)))
		  open(mindist(object))
		  for(j in seq_along(id)){
			  if(verbose) message("\t Pruning sample ", j, " of ", length(id))
			  sampleId <- id[j]
			  rd <- ranges[sampleNames(ranges) == sampleId, ]
			  stopifnot(nrow(rd) > 0)
			  ## assign the mad of the minimum distance to the ranges
			  k <- match(sampleId, sampleNames(object))
			  ##rd$mad <- object[[1]]$mindist.mad[k]
			  genomdat <- as.numeric(mindist(object)[, k])
			  rdList[[j]] <- pruneMD(genomdat,
						 rd,
						 physical.pos=position(object),  ##fD$position,
						 lambda=lambda,
						 MIN.CHANGE=min.change,
						 SCALE.EXP=scale.exp,
						 MIN.COVERAGE=min.coverage)
		  }
		  close(mindist(object))
		  if(length(rdList) == 1) {
			  rd <- rdList[[1]]
		  } else {
			  rdList <- rdList[!sapply(rdList, is.null)]
			  ##rdList <- lapply(rdList, function(x) as(x, "RangedData"))
			  rdl <- RangedDataList(rdList)
			  rd <- stack(rdl)
			  ##rd <- as(rd, "RangedDataCNV")
			  ix <- match("sample", colnames(rd))
			  if(length(ix) > 0) rd <- rd[, -ix]
		  }
		  return(rd)
	 })

## CIDR_Name is not a general label -- need to generalize these functions
setMethod("offspringNames", signature(object="TrioSet"), function(object){
	phenoData2(object)[, "CIDR_Name", "O"]
})
setMethod("fatherNames", signature(object="TrioSet"), function(object){
	phenoData2(object)[, "CIDR_Name", "F"]
})
setMethod("motherNames", signature(object="TrioSet"), function(object){
	phenoData2(object)[, "CIDR_Name", "M"]
})
fmoNames <- function(object){
	tmp <- cbind(fatherNames(object), motherNames(object), offspringNames(object))
	colnames(tmp) <- c("F", "M", "O")
	return(tmp)
}

logrPanel <- function(x, y, segments=TRUE,
		      range, fmonames,
		      cbs.segs,
		      md.segs, ..., subscripts){
	panel.grid(v=10,h=10, "grey")
	panel.xyplot(x, y, ...)
	if(segments){
		CHR <- range$chrom
		if(missing(cbs.segs)){
			loadRangesCbs(beadstudiodir(), pattern=paste("cbs_chr", CHR, sep=""), name="cbs.segs")
		}
		what <- switch(paste("p", panel.number(), sep=""),
			       p1="min dist",
			       p2="offspring",
			       p3="mother",
			       p4="father",
			       NULL)
		stopifnot(!is.null(what))
		if(what=="father")
			cbs.sub <- cbs.segs[cbs.segs$id==fmoNames[1], ]
		if(what=="mother")
			cbs.sub <- cbs.segs[cbs.segs$id==fmoNames[2], ]
		if(what=="offspring")
			cbs.sub <- cbs.segs[cbs.segs$id==fmoNames[3], ]
		if(what=="min dist"){
			if(!missing(md.segs)){
				cbs.sub <- md.segs[md.segs$id %in% ss(range$id), ]
				##cbs.sub <- dranges[substr(dranges$id, 1, 5) %in% id, ]
				cbs.sub$seg.mean <- -1*cbs.sub$seg.mean
			}
		}
		if(nrow(cbs.sub) > 0){
			cbs.sub$seg.mean[cbs.sub$seg.mean < ylimit[1]] <- ylimit[1] + 0.2
			cbs.sub$seg.mean[cbs.sub$seg.mean > ylimit[2]] <- ylimit[2] - 0.2
			stopifnot(nrow(cbs.sub) > 0)
			panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2,col="black")#gp=gpar("lwd"=2))
		}
	}
}

setMethod("xyplot", signature(x="formula", data="TrioSet"),
	  function(x, data, ...){
		  if(!"panel" %in% names(list(...))){
			  panel <- logrPanel
		  }
		  data <- Beaty:::todf(data, ...)
		  if("panelLabels" %in% names(list(...))){
			  panelLabels <- list(...)[["panelLabels"]]
			  data <- data[data$id %in% panelLabels, ]
		  }
		  xyplot(x=x, data=data,
			 panel=panel, ...)
			 ##layout=c(1,4),
			 ##index.cond=list(4:1),
			 ##pch=pch,
			 ##cex=cex,
			 ##xlab="Mb",
			 ##ylab="",
			 ##border="grey60",
			 ##range.object=rd[i,],
			 ##scales=list(x=list(tick.number=10, cex=cex.scale, tck=c(1,0)),
			 ##alternating=rep(1, 4),
			 ##y=list(cex=cex.scale, tck=c(1,0))),
			 ##par.strip.text=list(lines=0.8, cex=0.7),
			 ##alpha=alpha,
			 ##highlight=highlight,
			 ##col=col,
			 ##fill=fill, ...)##, main=rd$id[i]))
	  })
