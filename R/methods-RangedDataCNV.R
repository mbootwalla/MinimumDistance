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

setMethod("xprune", signature(object="TrioSetList", ranges="RangedDataCNV"),
	  function(object, ranges, ...){
		  rdList <- lapply(object, xprune, ranges=ranges, ...)
	  })

setMethod("xprune", signature(ranges="RangedDataCNV", trioSet="TrioSet"),
	  function(object, ranges, ...){
		  CHR <- unique(chromosome(object))
		  if(missing(id)) id <- unique(id(ranges))
		  index <- which(chromosome(ranges) == CHR & sampleNames(ranges) %in% id)
		  ranges <- ranges[index, ]
		  rdList <- vector("list", length(unique(id)))
		  open(mindist(object))
		  for(j in seq_along(id)){
			  sampleId <- id[j]
			  rd <- ranges[sampleNames(ranges) == sampleId, ]
			  stopifnot(nrow(rd) > 0)
			  ## assign the mad of the minimum distance to the ranges
			  k <- match(sampleId, sampleNames(object))
			  ##rd$mad <- object[[1]]$mindist.mad[k]
			  genomdat <- as.numeric(mindist(object)[, k])
			  rdList[[j]] <- prune(genomdat,
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
			  rd <- RangedDataList(rdList)
			  rd <- stack(rd)
			  ix <- match("sample", colnames(rd))
			  if(length(ix) > 0) rd <- rd[, -ix]
		  }
		  return(rd)
	 })

