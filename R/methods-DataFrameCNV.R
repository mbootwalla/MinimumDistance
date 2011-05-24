##setMethod("dimnames", "DataFrameCNV", function(x) list(x@row.names, names(x)))
##setMethod("dim", "DataFrameCNV", function(x) c(length(x@row.names), length(x)))

setMethod("xyplot", signature(x="formula", data="DataFrameCNV"),
	  function(x, data, ...){
		  if(!"panel" %in% names(list(...))){
			  panel <- my.xypanel
		  }
		  data(chromosomeAnnotation)
		  data <- as(data, "data.frame")
		  ## could we do UseMethod here?
		  fig <- xyplot(x, data=data,
				panel=panel,
				x0=data$x0,
				x1=data$x1,
				col=data$col,
				border=data$border,
				alpha=1,
				chr=data$chr,
				chr.size=data$chr.size,
				coverage=data$coverage,
				xlab="Mb",
				ylab="offspring index",
				##key=getKey(data),
				par.strip.text=list(lines=0.7, cex=0.6),
				prepanel=prepanel.fxn,
				max.y=max(data$y),
				chromosomeAnnotation=chromosomeAnnotation,...)
		  return(fig)
	  })

setMethod("[", signature(x="DataFrameCNV"),
	  function(x, i, j, ..., drop=FALSE){
		  xlist <- as(x, "data.frame")
		  xlist <- xlist[i, ]
		  x <- as(xlist, "DataFrameCNV")
		  return(x)
	  })

getKey <- function(df, space="top"){
##	labels <- factor(unique(df$statelabels), levels=c("loss (homozygous)",
##					      "loss (hemizygous)",
##					      "gain"), ordered=TRUE)
	labels <- unique(df$statelabels)
	labels <- labels[order(labels)]
	##col <- palette[as.integer(labels)]
	col <- unique(df$col)
	##col <- col[order(col)]
	##col <- df$col[match(states, df$state)]
	mykey <- simpleKey(as.character(labels), points=FALSE,
			   rectangles=TRUE, col=col, space=space)
	mykey$rectangles[["border"]] <- mykey$rectangles[["col"]] <- col
	mykey
}

my.xypanel <- function(x, y,
		       x0, x1, chr.size,
		       col, border, coverage,
		       chr, show.coverage=TRUE,
		       max.y,
		       chromosomeAnnotation,
		       addCentromere=TRUE,
		       ..., subscripts){
	panel.grid(h=-1, v=10)
	panel.xyplot(x, y, ..., subscripts)
	h <- 0.75
	lrect(xleft=x0[subscripts],
	      xright=x1[subscripts],
	      ybottom=y-h/2,
	      ytop=y+h/2,
	      border=border[subscripts],
	      col=col[subscripts], ...)
	if(show.coverage)
		ltext(x, y,labels=coverage[subscripts], cex=0.6)
	##plot centromere
	if(addCentromere){
		chr <- unique(as.integer(as.character(chr)))
		coords <- chromosomeAnnotation[chr, 1:2]/1e6
		lrect(xleft=coords[1],
		      xright=coords[2],
		      ybottom=0,
		      ytop=max.y+h/2,
		      col="grey",
		      border="grey")
	}
}

prepanel.fxn <- function(x,y, chr.size, ..., subscripts){
	list(xlim=c(0, unique(chr.size[subscripts])), ylim=range(as.integer(as.factor(y[subscripts]))))
}
