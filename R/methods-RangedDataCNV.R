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

setMethod("xprune", signature(ranges="RangedDataCNV"),
	  function(ranges, id, trioSets, lambda=0.05,
		   min.change=0.02, min.coverage=3,
		   scale.exp, CHR,
		   verbose, ...){
		 trioSet <- trioSets[[CHR]]
		 if(missing(id)) id <- unique(id(ranges))
		 index <- which(chromosome(ranges) == CHR & sampleNames(ranges) %in% id)
		 ranges <- ranges[index, ]
		 rdList <- vector("list", length(unique(id)))
		 open(mindist(trioSet))
		 for(j in seq_along(id)){
			 sampleId <- id[j]
			 rd <- ranges[sampleNames(ranges) == sampleId, ]
			 stopifnot(nrow(rd) > 0)
			 ## assign the mad of the minimum distance to the ranges
			 k <- match(sampleId, sampleNames(trioSets))
			 rd$mad <- trioSets[[1]]$mindist.mad[k]
			 ##chrom <- unique(chromosome(rd))
			 ##for(l in seq_along(chrom)){
				## CHR <- chrom[l]
			 ##rd2 <- rd[chromosome(rd) == CHR, ]
			 ##stopifnot(nrow(rd2) > 0)
			 ##trioSet <- trioSets[[CHR]]
			 genomdat <- as.numeric(mindist(trioSet)[, k])
			 rdList[[j]] <- prune(genomdat,
					      rd,
					      physical.pos=position(trioSet),  ##fD$position,
					      lambda=lambda,
					      MIN.CHANGE=min.change,
					      SCALE.EXP=scale.exp,
					      MIN.COVERAGE=min.coverage)
		 }
		 if(length(rdList) == 1) {
			 rd <- rdList[[1]]
		 } else {
			 rdList <- rdList[!sapply(rdList, is.null)]
			 rd <- RangedDataList(rdList)
			 rd <- stack(rd)
			 ##rd <- do.call("c", rdList)
		 }
##		 range.pruned <- RangedDataCNV(IRanges(start(rd), end(rd)),
##					       chrom=rd$chrom,
##					       id=rd$id,
##					       num.mark=rd$num.mark,
##					       seg.mean=rd$seg.mean,
##					       start.index=rd$start.index,
##					       end.index=rd$end.index,
##					       mindist.mad=rd$mad)
	 })

