setMethod("state", signature(object="RangedDataCNV"), function(object) object$state)
setMethod("nMarkers", signature(object="RangedDataCNV"), function(object) object$num.mark)
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

setMethod("plot", signature(x="RangedDataCNV", y="missing"),
	  function(x, y, ...){
		  df <- todf(x)
		  plot(x=df, ...)
	  })

setMethod("todf", signature(object="RangedDataCNV"), function(object, ...){
	require(SNPchip)
	data(chromosomeAnnotation)
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
			  col=as.integer(as.factor(state(object))),
			  state=state(object),
			  y=y,
			  stringsAsFactors=FALSE)
	dat$chr <- as.factor(dat$chr)
	dat <- new("DataFrameCNV", dat)
	return(dat)
})


setMethod("plot", signature(x="DataFrameCNV"),
	  function(x, ...){
		  data(chromosomeAnnotation)
		  df <- as.data.frame(object@.Data)
		  colnames(df) <- names(object)
		  df$x <- df$midpoint
		  fig <- xyplot(y~x, data=df,
				panel=my.xypanel,
				x0=df$x0,
				x1=df$x1,
				col=df$col,
				border=df$border,
				alpha=1,
				chr=df$chr,
				chr.size=df$chr.size,
				coverage=df$coverage,
				xlab="Mb",
				ylab="offspring index",
				key=getKey(df),
				par.strip.text=list(lines=0.7, cex=0.6),
				prepanel=prepanel.fxn,
				max.y=max(df$y),
				chromosomeAnnotation=chromosomeAnnotation,
				...)
		  return(fig)
	  })

##setMethod("rbind", "RangedDataCNV", function(..., deparse.level=1){
##	x <- lapply(list(...), function(x) as(x, "RangedData"))
##	rdl <- RangedDataList(x)
##	rd <- IRanges:::stack(rdl)
##	##rd <- rd[, -ncol(rd)]
##	rd <- as(rd, "RangedDataCNV")
##	rd
##})



