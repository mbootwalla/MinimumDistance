RangedDataCNV <- function(ranges=IRanges(),
			  values,
			  start,
			  end,
			  chromosome,
			  coverage,
			  sampleId,
			  startIndexInChromosome,
			  endIndexInChromosome,
			  ...){
	if(!missing(ranges) & !missing(values)){
		object <- new("RangedDataCNV", ranges=ranges, values=values)
		return(object)
	}
	if(!missing(end) && !missing(start))
		ranges <- IRanges(start, end)
	if(missing(chromosome))
		chromosome <- vector("integer", length(ranges))
	if(missing(coverage))
		coverage <- vector("integer", length(ranges))
	if(missing(sampleId))
		sampleId <- vector("character", length(ranges))
	if(missing(startIndexInChromosome))
		startIndexInChromosome <- vector("integer", length(ranges))
	if(missing(endIndexInChromosome))
		endIndexInChromosome <- vector("integer", length(ranges))
	rd <- RangedData(ranges,
			 chrom=chromosome,
			 num.mark=coverage,
			 id=sampleId,
			 start.index=startIndexInChromosome,
			 end.index=endIndexInChromosome, ...)##, ...)
	new("RangedDataCNV", ranges=ranges(rd), values=IRanges:::values(rd))
}

RangedDataCBS <- function(ranges=IRanges(),
			  seg.mean=vector("numeric", length(ranges)), ...){
	rd <- RangedDataCNV(ranges=ranges, seg.mean=seg.mean, ...)
	new("RangedDataCBS", ranges=ranges(rd), values=values(rd))
}
RangedDataCBS2 <- function(ranges=IRanges(),
			   state=vector("character", length(ranges)), ...){
	rd <- RangedDataCBS(ranges=ranges, state=state, ...)
	new("RangedDataCBS2", ranges=ranges(rd), values=values(rd))
}


RangedDataHMM <- function(ranges=IRanges(),
			  state=vector("integer", length(ranges)), ...){
	rd <- RangedDataCNV(ranges=ranges, state=state, ...)
	new("RangedDataHMM", ranges=ranges(rd), values=values(rd))
}

setMethod("state", signature(object="RangedDataCNV"), function(object) object$state)
setMethod("nMarkers", signature(object="RangedDataCNV"), function(object) object$num.mark)
setMethod("coverage", signature(object="RangedDataCNV"), function(object) object$num.mark)
setMethod("sampleNames", signature(object="RangedDataCNV"), function(object) object$id)
##setMethod("RangedDataCNV", signature(ranges="IRanges"),
##	  function(ranges=IRanges(), ...,
##		   space=NULL,
##		   universe=NULL){
##		  nms <- names(list(...))
##		  stopifnot(c("chrom", "id", "num.mark", "seg.mean", "start.index", "end.index") %in% nms)
##		  rd <- RangedData(ranges=ranges, ..., space=space, universe=universe)
##		  rd2 <- as(rd, "RangedDataCNV")
##		  return(rd2)
##	  })

setAs("RangedData", "RangedDataCBS", function(from){
	RangedDataCBS(ranges=ranges(from),
		      values=values(from))
})
setAs("RangedData", "RangedDataHMM", function(from){
	RangedDataHMM(ranges=ranges(from),
		      values=values(from))
})

setMethod("chromosome", signature(object="RangedDataCNV"), function(object) object$chrom)

##setMethod("plot", signature(x="RangedDataCNV", y="missing"),
##	  function(x, y, ...){
##		  df <- todf(x)
##		  plot(x=df, ...)
##	  })

setMethod("todf", signature(object="RangedDataCNV"), function(object, col=1:3, verbose=TRUE, ordered.y, ...){
	require(SNPchip)
	data(chromosomeAnnotation)
	if(!"col" %in% names(list(...))) col <- 1:3 else col <- list(...)[["col"]]
	##object <- object[object$state != normalIndex(hmm.params), ]
	is.denovo <- isDenovo(state(object))
	if(verbose) message("Dropping ranges that are not denovo")
	object <- object[is.denovo, ]
	h <- 0.75
	meanSegment <- apply(cbind(start(object), end(object)), 1, mean)
	chr.size <- chromosomeAnnotation[1:22, "chromosomeSize"]
	chrom <- chromosome(object)
	chr.size <- chr.size[chrom]
	if(missing(ordered.y)){
		y <- split(sampleNames(object), chrom)
		y <- lapply(y, function(x){
			tmp <- as.numeric(as.factor(x))
			names(tmp) <- as.character(x)
			tmp
		})
		y <- unlist(y)
		nms2 <- paste(chromosome(object), sampleNames(object), sep=".")
		if(!identical(names(y), nms2)){
			y <- y[match(nms2, names(y))]
			stopifnot(identical(names(y), nms2))
		}
	} else {
		y <- ordered.y
	}
	##
	states <- state(object)
	xstates <- rep(NA, length(states))
	xstates[states %in% offspring.homozygous()] <- "loss (homozygous)"
	xstates[states %in% offspring.hemizygous()] <- "loss (hemizygous)"
	xstates[states %in% duplicationStates()] <- "gain"
	##
	labels <- factor(xstates, levels=c("loss (homozygous)",
				  "loss (hemizygous)",
				  "gain"), ordered=TRUE)
	colors <- col[as.integer(as.factor(labels))]
	##colors <- factor(colors, levels=col, ordered=TRUE)
	##
	dat <- data.frame(x0=start(object)/1e6,
			  x1=end(object)/1e6,
			  y0=y-h/2,
			  y1=y+h/2,
			  chr=chromosome(object),
			  coverage=nMarkers(object),
			  midpoint=meanSegment/1e6,
			  id=sampleNames(object),
			  chr.size=chr.size/1e6,
			  border=rep("black", nrow(object)),
			  col=colors,
			  statelabels=labels,
			  state=state(object),
			  y=y,
			  stringsAsFactors=FALSE)
	dat$chr <- as.factor(dat$chr)
	dat <- as(dat, "DataFrameCNV")
	##dat <- new("DataFrameCNV", dat)
	return(dat)
})

setMethod("rbind", "RangedDataCNV",
	  function(..., deparse.level=1){
		  callNextMethod(..., deparse.level=deparse.level)
	  })


##setMethod("[", signature(x="RangedDataCNV"),
##	  function(x, i, j, ..., drop=FALSE){
##		  ## The "[" method for RangedData does not have separate methods for the slot elements
##		  ## It seems that the easiest way to subset an object extending the RangedData class is to use
##		  ## the RangeData method directly -- this will return an object of class RangedData
##		  rd <- callNextMethod(x, i, j, ..., drop=drop)
##		  ## Now, update the components of 'x'
##		  x@ranges <- ranges(rd)
##		  x@values <- values(rd)
##		  return(x)
##	  })



