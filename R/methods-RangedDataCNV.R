setMethod("state", signature(object="RangedDataCNV"), function(object) object$state)
setMethod("nMarkers", signature(object="RangedDataCNV"), function(object) object$num.mark)
setMethod("coverage", signature(object="RangedDataCNV"), function(object) object$num.mark)
setMethod("sampleNames", signature(object="RangedDataCNV"), function(object) object$id)
setMethod("RangedDataCNV", signature(ranges="IRanges"),
	  function(ranges=IRanges(), ...,
		   space=NULL,
		   universe=NULL){
		  nms <- names(list(...))
		  stopifnot(c("chrom", "id", "num.mark", "seg.mean", "start.index", "end.index") %in% nms)
		  rd <- RangedData(ranges=ranges, ..., space=space, universe=universe)
		  rd2 <- as(rd, "RangedDataCNV")
		  return(rd2)
	  })

setMethod("chromosome", signature(object="RangedDataCNV"), function(object) object$chrom)

##setMethod("plot", signature(x="RangedDataCNV", y="missing"),
##	  function(x, y, ...){
##		  df <- todf(x)
##		  plot(x=df, ...)
##	  })




setMethod("todf", signature(object="RangedDataCNV"), function(object, col=1:3, ...){
	require(SNPchip)
	data(chromosomeAnnotation)
	if(!"col" %in% names(list(...))) col <- 1:3 else col <- list(...)[["col"]]
	##object <- object[object$state != normalIndex(hmm.params), ]
	is.denovo <- isDenovo(state(object))
	object <- object[is.denovo, ]
	h <- 0.75
	meanSegment <- apply(cbind(start(object), end(object)), 1, mean)
	chr.size <- chromosomeAnnotation[1:22, "chromosomeSize"]
	chrom <- chromosome(object)
	chr.size <- chr.size[chrom]
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

##setMethod("rbind", "RangedDataCNV", function(..., deparse.level=1){
##	x <- lapply(list(...), function(x) as(x, "RangedData"))
##	rdl <- RangedDataList(x)
##	rd <- IRanges:::stack(rdl)
##	##rd <- rd[, -ncol(rd)]
##	rd <- as(rd, "RangedDataCNV")
##	rd
##})



