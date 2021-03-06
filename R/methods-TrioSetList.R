setMethod("xsegment", signature(object="TrioSetList"),
	  function(object, id, ..., verbose=FALSE, DNAcopy.verbose=0){
		  dfl <- lapply(object, xsegment, id=id, ..., verbose=verbose,DNAcopy.verbose=DNAcopy.verbose)
		  df <- do.call("rbind", dfl)
		  return(df)
	  })
setMethod("sampleNames", signature(object="TrioSetList"),
	  function(object) sampleNames(object[[1]]))
setMethod("ncol", signature(x="TrioSetList"),
	  function(x) ncol(x[[1]]))
setMethod("nrow", signature(x="TrioSetList"),
	  function(x) nrow(x[[1]]))
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
		  return(rdList)
	  })



setMethod("offspringNames", signature(object="TrioSetList"), function(object) offspringNames(object[[1]]))
setMethod("fatherNames", signature(object="TrioSetList"), function(object) fatherNames(object[[1]]))
setMethod("motherNames", signature(object="TrioSetList"), function(object) motherNames(object[[1]]))

setMethod("computeBayesFactor", signature(object="TrioSetList"),
	  function(object, ranges, id, states, baf.sds, mu.logr,
		   log.pi, tau, normal.index, a,
		   prOutlier=c(0.01, 1e-5),
		   prMosaic=0.01,
		   prob.nonMendelian,
		   df0,
		   verbose,
		   returnEmission){
		  if(missing(id)) id <- unique(ranges$id) else stopifnot(id %in% unique(ranges$id))
		  chromosomes <- sapply(object, function(x) unique(chromosome(x)))
		  ranges <- ranges[chromosome(ranges) %in% chromosomes, ]
		  ranges <- ranges[ranges$id %in% id, ]
##		  if(!"bayes.factor" %in% colnames(ranges)){
##			  ranges$bayes.factor <- NA
##		  }
		  if(!"lik.state" %in% colnames(ranges)){
			  ranges$lik.state <- NA
		  }
		  if(!"lik.norm" %in% colnames(ranges)){
		  	  ranges$lik.norm <- NA
		  }
		  if(!"argmax" %in% colnames(ranges)){
			  ranges$argmax <- NA
		  }
		  for(i in seq_along(object)){
			  if(verbose)
				  message("\tProcessing chromosome ", i, " of ", length(object))
			  CHR <- unique(chromosome(object[[i]]))
			  j <- which(chromosome(ranges) == CHR)
			  if(length(j) < 1) next()
			  rd <- computeBayesFactor(object[[i]],
						   ranges[j, ],
						   id=id,
						   states=states,
						   baf.sds=baf.sds,
						   mu.logr=mu.logr,
						   log.pi=log.pi,
						   tau=tau,
						   normal.index=normal.index,
						   a=a,
						   prOutlier=prOutlier,
						   prMosaic=prMosaic,
						   prob.nonMendelian=prob.nonMendelian,
						   df0=df0,
						   returnEmission=returnEmission,
						   verbose=verbose)
			  if(returnEmission) return(rd)
			  ranges$lik.state[j] <- rd$lik.state
			  ranges$argmax[j] <- rd$argmax
			  ranges$lik.norm[j] <- rd$lik.norm
			  ##ranges$DN[j] <- rd$DN
		  }
		  return(ranges)
	  })

setMethod("[", signature(x="TrioSetList"),
	  function(x, i, j, ..., drop=FALSE){
		  xlist <- as(x, "list")
		  xlist <- xlist[i]
		  x <- as(xlist, "TrioSetList")
		  return(x)
	  })

setMethod("xyplot", signature(x="formula", data="TrioSetList"),
	  function(x, data, ...){
		  stopifnot("range" %in% names(list(...)))
		  range <- list(...)[["range"]]
		  stopifnot(nrow(range)==1)
		  trioSet <- data[[range$chrom]]
		  xyplot(x, trioSet, ...)
	  })
