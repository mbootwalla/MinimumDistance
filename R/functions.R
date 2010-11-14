scaleSnr <- function(snr, min.scale, max.scale){
	tmp <- (snr-min(snr))/(max(snr)-min(snr))  ## 0 -> 1
	b <- 1/(max.scale-min.scale)
	a <- min.scale*b
	bg.scale <- (tmp + a)/b
	return(bg.scale)
}

isDeletion <- function(x){
	if(length(grep("-", x)) > 0){
		tmp <- strsplit(x, "_")[[1]]
		state <- substr(tmp, 3, 3)
		state <- ifelse(any(state < 3), TRUE, FALSE)
	} else{
		state <- as.integer(substr(x, 3, 3))
		state <- ifelse(state < 3, TRUE, FALSE)
	}
	state
}

readBiparentalMatrix <- function(tmpdir="/tmp/bpi"){
	biparental <- read.csv(file.path(tmpdir, "numberInformative.csv"), as.is=TRUE)
	##convert columns 2-4 to an integer
	for(j in 2:4) biparental[, j] <- suppressWarnings(as.integer(biparental[, j]))
	bpi <- matrix(unlist(biparental[, 2:4]), nrow(biparental), 3)
	exclude <- which(rowSums(is.na(bpi)) == 3)
	fns <- biparental[-exclude, 1]
	bpi <- bpi[-exclude, ]
	colnames(bpi) <- colnames(biparental)[2:4]
	rownames(bpi) <- fns
	delta <- bpi[, "nBiparental"] - bpi[, "nNotBiparental"]
	bpi <- cbind(bpi, delta)
	colnames(bpi)[4] <- "delta"
	bpi
}

constructDisjointRangedData <- function(rD, denovoSet, ids.exclude, minoverlap=25e3){
	numberMarkers <- disjointRanges <- nshared <- vector("list", 22)
	for(CHR in 1:22){
		cat(CHR, " ")
		i <- which(rD$chrom==CHR & !(rD$id %in% ids.exclude))
		subject <- IRanges(start(rD)[i], end(rD)[i])
		subject.tree <- IntervalTree(subject)
		##construct query as the disjoint ranges
		starts <- sort(union(start(subject.tree), end(subject.tree)))
		disjointRanges[[CHR]] <- IRanges(starts[-length(starts)], starts[-1])
		chrom <- chromosome(denovoSet)
		pos <- position(denovoSet)
		pos <- pos[chrom == CHR & !is.na(chrom)]
		query.markers <- IRanges(start=pos-12,
					 end=pos+12)
		tree.markers <- IntervalTree(query.markers)
		disjoint.tree <- IntervalTree(disjointRanges[[CHR]])

		##how many markers are in each segment
		numberMarkers[[CHR]] <- countOverlaps(disjoint.tree, tree.markers)

		##specify a minimum size of overlap (e.g., 50kb)
		nshared[[CHR]] <- countOverlaps(disjointRanges[[CHR]], subject.tree, minoverlap=minoverlap)
		if(CHR==22) cat("\n")
	}

	##the nshared corresponds to the number of denovo regions in disjoint bin
	##once we find the peak, we could define a region as the union of all
	##segments in which there is at least one overlap.
	tmp <- do.call("c", disjointRanges)
	chrom <- rep(1:22, sapply(nshared, length))
	nshared <- unlist(nshared)
	numberMarkers <- unlist(numberMarkers)
	rangedData <- RangedData(tmp, numberMarkers=numberMarkers, numberOverlaps=nshared,
				 chrom=chrom)
	rangedData <- rangedData[order(rangedData$numberOverlaps, rangedData$numberMarkers, decreasing=TRUE), ]
	return(rangedData)
}

summarizeBias <- function(bias2, cnSet){
	bias2.normal <- bias2[, cnSet$trisomy==0, ]
	bias2.trisomy <- bias2[, cnSet$trisomy==1, ]
	bias2.snps <- cbind(rowMedians(bias2.normal[, , "snps"], na.rm=TRUE),
			    rowMedians(bias2.trisomy[, , "snps"], na.rm=TRUE))
	bias2.nps <- cbind(rowMedians(bias2.normal[, , "nps"], na.rm=TRUE),
			   rowMedians(bias2.trisomy[, , "nps"], na.rm=TRUE))
	colnames(bias2.snps) <- c("CN2", "CN3")
	colnames(bias2.nps) <- c("CN2", "CN3")
	bias2.overall <- bias2.snps + bias2.nps
	colnames(bias2.overall) <- c("CN2", "CN3")
	bias2.overall[1, ] <- bias2.snps[1, ]
	bias2.marginal <- bias2.overall[, 1]+bias2.overall[, 2]
	bias2.marginal[1] <- bias2.overall[1,1]
	return(list(marginal=bias2.marginal,
		    overall=bias2.overall))
}
callDenovoEvents <- function(denovoSet, bsSet, trios.index, THR){
	segmentMeans <- assayData(bsSet)[["segmentMeans"]]
	trios.index <- getTriosIndex(bsSet)
	open(segmentMeans)
        D <- assayData(denovoSet)[["exprs"]]
        open(D)
        k <- 1
##	logR <- assayData(bsSet)[["logR"]]
##	open(logR)
##	mads <- vector("numeric", ncol(bsSet))
##	ii <- which(chromosome(bsSet) < 23)
##	for(j in 1:ncol(bsSet)) mads[j] <- mad(logR[ii, j], na.rm=TRUE)
##	mad.cutoff <- quantile(mads, probs=0.99)
	for(j in 1:ncol(denovoSet)){
		offspring.id <- sampleNames(denovoSet)[j]
		family.index <- which(bsSet$family == denovoSet$family[j])
		sns <- bsSet$Sample.Name[family.index]
		individ <- paste("0", bsSet$individ[family.index], sep="")
		whoisit <- sapply(individ, who)
		condition <- any("offspring" %in% whoisit) & any("father" %in% whoisit) & any("mother" %in% whoisit)
		if(!condition){
			D[, j] <- 0
			next()
		}
                if(j %% 10 == 0) cat(j, " ")
		parent.index <- c(family.index[whoisit == "mother"], family.index[whoisit == "father"])
		if(length(parent.index) < 2) {
			message("parent.index < 2 for " , offspring.id)
			D[, j] <- 0
			next()
		}
		offspring.index <- family.index[whoisit == "offspring"]

		##check that both parents have normal copy number
		segs.parents <- segmentMeans[, parent.index, drop=FALSE]
                parentsBothNormal <- rowSums(segs.parents > THR[[2]], na.rm=TRUE) == 2
                ##denovo indicator
                offspringHasDeletion <- segmentMeans[, offspring.index] < THR[[1]]
                D[, j] <- offspringHasDeletion & parentsBothNormal
		## Check the variance of the log R values
##		index <- which(D[, j] == 1)
##		mad.in.region <- apply(logR[index, parent.index], 2, mad, na.rm=TRUE)
##		##parent variance should be small
##		D[index, j] <- ifelse(mad.in.region
        }
        close(D)
        close(segmentMeans)
        return(TRUE)
}


getIndex <- function(query, bsSet, window.size=1e6){
	ix <- match(query$id, sampleNames(bsSet))
	family.index <- which(bsSet$family == bsSet$family[ix])
	row.index <- position(bsSet) >= (start(query) - window.size) & position(bsSet) <= end(query) + window.size
	row.index <- row.index & chromosome(bsSet) == query$chrom
	row.index <- which(row.index)

	sns <- sampleNames(bsSet)[family.index]
	individ <- paste("0", bsSet$individ[family.index], sep="")
	whoisit <- sapply(individ, who)
	family.index <- family.index[order(whoisit)]
	list(row.index, family.index)
}







plotImage <- function(denovoSet, indices, query, minoverlap){
	require(SNPchip)
	chrom <- query$chrom
        pathto <- system.file("hg18", package = "SNPchip")
        cytoband <- read.table(file.path(pathto, "cytoBand.txt"),
			       as.is = TRUE)
        colnames(cytoband) <- c("chrom", "start", "end", "name",
				"gieStain")
	data(chromosomeAnnotation)
	marker.index <- indices[[1]]
	col.index <- indices[[2]]
	M <- copyNumber(denovoSet)[marker.index, col.index]
	x <- position(denovoSet)[marker.index]
	xlim <- range(x)
	##image of region
	col.image <- c("white", "black")
	centromere.coords <- chromosomeAnnotation[query$chrom, ]
	layout(mat=matrix(c(1, 1,
	       2, 3,
	       4, 4), ncol=2, byrow=TRUE), heights=c(0.1, 0.8, 0.1), widths=c(0.02, 0.95))
	hist(x, freq=TRUE, breaks=length(x)/2, xlim=xlim, xaxs="i", main="", xaxt="n")
	rug(x, side=3, col="light blue")
	sns <- substr(colnames(M), 1, 8)
	ylim=c(-5, ncol(M)+5)
	label.cols <- grey(seq(0, by=1/ncol(denovoSet)))
	label.cols <- label.cols[col.index]
	label.matrix <- matrix(col.index, nrow=1, byrow=FALSE)
	image(0, 1:ncol(M), z=label.matrix, col=label.cols, ylim=ylim, ylab="", yaxt="n", xaxt="n")
	image(x, 1:ncol(M), M, col=col.image, ylab="", xlim=xlim, xaxs="i",
	      xlab="Mb", yaxt="n", bg="lightblue",
	      xaxt="n", ylim=ylim)
	abline(v=c(start(query), end(query)), col="blue", lty=2)
	##axis(3, at=x, labels=FALSE)
	## how to draw a second vertical bar that shows the union?
	if(ncol(M) > 50){
		axis(2, at=pretty(1:ncol(M)), labels=pretty(1:ncol(M)), cex.axis=0.8)
	} else {
		axis(2, at=1:ncol(M), labels=sns, cex.axis=0.6)
	}
	axis(1, at=pretty(xlim), labels=pretty(xlim)/1e6)
	## show what minoverlap looks like
	xx <- c(mean(x), mean(x)+minoverlap)
	graphics:::segments(x0=xx[1], x1=xx[2], y0=-3, y1=-3, col="blue")
	graphics:::segments(x0=xx[1], x1=xx[1], y0=-4, y1=-2, col="blue")
	graphics:::segments(x0=xx[2], x1=xx[2], y0=-4, y1=-2, col="blue")
	text(x=xx[2]+minoverlap, y=-3, labels=paste(minoverlap/1e3, "kb"), cex=0.7, adj=0)
	invisible(plotCytoband(query$chrom, label.cytoband=FALSE, cytoband.ycoords=c(0.75, 1.25)))
	graphics:::segments(x0=xlim[1], x1=xlim[2], y0=0.5, y1=0.5, col="blue")
	graphics:::segments(x0=xlim[1], x1=xlim[1], y0=0.45, y1=0.55, col="blue")
	graphics:::segments(x0=xlim[2], x1=xlim[2], y0=0.45, y1=0.55, col="blue")
	mtext(paste("chr", query$chrom), cex=1.5, side=1)
}

findOffspringParentIndex <- function(indices, bsSet, denovoSet){
	sampleIndex <- indices[[2]]
	bsSetIndex <- which(bsSet$family %in% denovoSet$family[sampleIndex])
	bsSetIndexByFamily <- split(bsSetIndex, bsSet$family[bsSetIndex])
	bsSetIndexByFamily
}


plotTrioLogRInRegion <- function(marker.index,
				 sample.index,
				 query, bsSet,
				 denovoSet,
				 window.size,...){
	lR <- logR(bsSet)[marker.index, sample.index]
	x <- position(bsSet)[marker.index]
	for(j in seq_along(sample.index)){
		plot(x, lR[, j], ...)
	}
}

callDenovoRegionsFromLogR <- function(object){
	 ##---------------------------------------------------------------------------
	 ##
	 ##  initialize object for storing denovo indicator
	 ##
	 ##---------------------------------------------------------------------------
	 if(!file.exists(file.path(outdir, "denovoSet2.rda"))){
		 denovoSet <- initializeEmptyDenovoSet(bsSet)
		 save(denovoSet, file=file.path(outdir, "denovoSet2.rda"))
	 } else {
		 load(file.path(outdir, "denovoSet2.rda"))
	 }
	 ##---------------------------------------------------------------------------
	 ##
	 ##  Call denovo events
	 ##
	 ##---------------------------------------------------------------------------
	 callDenovoEvents(denovoSet,
			  bsSet,
			  THR=c(log(1.25/2), log(1.75/2)))
	 ##---------------------------------------------------------------------------
	 ##
	 ##  Create RangedData objects from the denovo indicator matrix
	 ##
	 ##---------------------------------------------------------------------------
	 ocProbesets(20e3)
	 D <- assayData(denovoSet)[["exprs"]]
	 open(D)
         rdList <- vector("list", ncol(denovoSet))
         for(j in 1:ncol(denovoSet)){
                 cat(j, " ")
                 rdList[[j]] <- binary2RangedData(x=D[, j],
						  pos=position(denovoSet),
						  id=sampleNames(denovoSet)[j],
						  chr=chromosome(denovoSet))
         }
	 autosome.subset <- function(object) object <- object[object$chrom <= 22, ]
	 tmpList <- rdList[sapply(rdList, nrow) > 0]
	 rdList <- lapply(tmpList, autosome.subset)

         save(rdList, file=file.path(outdir, "denovo_rdList.rda")) ## save this since sometimes the next step doesn't work
         rD <- do.call("c", rdList)
         save(rD, file=file.path(outdir, "denovo_RangedData.rda"))
	 NULL
 }




checkRule <- function(query, ratioSet, sample.index){
         chr <- query$chrom
         start <- start(query)
         end <- end(query)
	 if("id" %in% colnames(query)){
		 id <- query$id
	 } else {
		 ## We need to find the ids for the offspring that contribute to the region being called denovo
		 if(missing(sample.index)) stop("if id is not in the query object, the sample.index must be specified")
		 if(length(sample.index) > 1){
			 id <- sample(names(sample.index), 1)
		 } else id <- sample.indexy

	 }
         pos <- position(ratioSet)
         chrom <- chromosome(ratioSet)
         family <- pData(ratioSet)[match(id, ratioSet$Sample.Name), "family"]
         family.index <- grep(family, pData(ratioSet)$family)
         individ <- paste("0", ratioSet$individ[family.index], sep="")
	 familynames <- sapply(individ, who)
         start.index <- min(which(pos == start & chrom == chr))
         end.index <- max(which(pos == end & chrom == chr))
         j <- start.index-50
         J <- end.index+50
	 rows <- (j:J)[chrom[j:J]==chr]
	 if(class(ratioSet)== "MultiSet"){
		 y <- assayData(ratioSet)[["logR"]][rows, family.index]
	 } else{
		 y <- logR(ratioSet)[rows, family.index]
	 }
	 ylim <- c(-2,1)
	 y[y < min(ylim)] <- min(ylim)
	 y[y > max(ylim)] <- max(ylim)
	 x <- pos[rows]
	 inRegion <- x >= start & x <= end
	 befRegion <- x < start
	 aftRegion <- x > end
         par(mfrow=c(3,1), las=1)
         for(k in seq_along(family.index)){
                 plot(x, y[, k], pch=21, cex=0.7, ylim=c(-2, 1))
		 legend("topright", bty="n", legend=familynames[k])
		 abline(v=c(start, end), lty=2, col="blue")
		 graphics:::segments(start, median(y[inRegion, k]), end, median(y[inRegion, k]), col="blue", lwd=2)
		 abline(h=0, col="grey70")
         }
 }





plotNumberDenovo <- function(nevents, nbases){
	## Drop samples that have an unusually large number of denovo events
	##pdf("suspicious.pdf", width=8, height=7)
	par(las=1, mfrow=c(2,2), mar=rep(0.5, 4), oma=c(4,4,2, 2))
	hist(log10(nevents), breaks=100, xaxt="n", main="")
	abline(v=2.5, col="blue")
	hist(log10(nbases), breaks=100, main="")
	abline(v=7, col="blue")
	mtext("log10(nbases)", line=3)

	I <- log10(nevents) >= 2.5 | log10(nbases) > 7
	plot(log10(nevents), log10(nbases), pch=21, cex=0.8)
	points(log10(nevents[I]),
	       log10(nbases[I]),
	       pch=21, bg="green3", cex=0.8)
	mtext("log10(nevents)", 1, line=3)
}

setAs("MultiSet", "LogRatioSet",
      function(from, to){
	      phenodata=phenoData(from)
	      sampleNames(phenodata)=phenodata$Sample.Name
	      prD <- protocolData(from)
	      sampleNames(prD) <- sampleNames(phenodata)
	      lR <- assayData(from)[["logR"]]
	      ix <- match(colnames(lR), sampleNames(phenodata))
	      phenodata <- phenodata[ix, ]
	      prD <- prD[ix, ]
	      new("LogRatioSet",
		  logRRatio=assayData(from)[["logR"]],
		  phenoData=phenodata,
		  protocolData=prD,
		  featureData=featureData(bsSet),
		  annotation=annotation(bsSet))
      })

subsetRangedData <- function(object, minimumOverlaps=10, excludeCentromeric=TRUE){
	rangedData <- rangedData[which(rangedData$numberOverlaps >= minimumOverlaps), ]
	if(excludeCentromeric){
		rangedData$overlapsCentromere <- rep(NA, nrow(rangedData))
		for(chr in unique(rangedData$chrom)){
			index <- which(rangedData$chrom == chr)
			tmp <- rangedData[index, ]
			rangeddata.tree <- IntervalTree(IRanges(start(tmp),end(tmp)))
			overlaps <- countOverlaps(rangeddata.tree,
						  chrAnn[chr, ])
			rangedData$overlapsCentromere[index] <- overlaps > 0
		}
		rangedData <- rangedData[which(!rangedData$overlapsCentromere), ]
	}
	return(rangedData)
}

getLrSet <- function() {
	load(file.path(outdir, "bsSet.rda"))
	bsSet <- get("bsSet")
	as(bsSet, "LogRatioSet")
}

doSegmentation <- function(){
	if(!exists("batch")) stop("batch variable should be defined in the submit_logR_cbs.R script")
	if(!exists("NN")) stop("batch size variable (NN) should be specified in submitter script")
	sample.index <- splitIndicesByLength(1:7599, NN)[[batch]]
	rD <- cbs(lrSet, sample.index=sample.index)
	save(rD, file=file.path(outdir, paste("rD_", batch, ".rda", sep="")))
	q("no")
}

guessChrom <- function(object){
	warning("Assumes segmentation was done in the order chromosome 1, chromosome 2, ...")
	starts <- start(object)
	chrom <- c(0, cumsum(diff(start(object)) < 0))+1
	if(length(unique(chrom)) > 25)
		stop("Probably did not work -- too many negative differences in the start coordinates")
	return(chrom)
}

getSegMeanRanges <- function(outdir){
	rd.fns <- list.files(file.path(outdir), pattern="rD_", full.names=TRUE)
	rdList <- vector("list", length(rd.fns))
	for(i in seq_along(rd.fns)){
		cat(i, " ")
		load(rd.fns[i])
		rD <- get("rD") ## list.  each element is a sample
		tmp <- rD[sapply(rD, class) == "RangedData"]
		tmp <- do.call("c", tmp)  ##create single RangedData object
##		if(!"chrom" %in% colnames(tmp)){
##			print("Chromosome not available in early version of cbs...taking a guess")
##			uid <- unique(tmp$id)
##			tmp$chrom <- NA
##			for(j in seq_along(uid)){
##				index <- which(tmp$id == uid[j])
##				suppressWarnings(tmp$chrom[index] <- guessChrom(tmp[index, ]))
##			}
##		}
		rdList[[i]] <- tmp
		if(i == length(rd.fns)) cat("\n")
	}
	rm(tmp); gc()
	##save(rdList, file=file.path(outdir, "rdList.rda"))
	tmp <- do.call("c", rdList)
	sns <- as.factor(tmp$id)
	seg.mean <- as.integer(tmp$seg.mean*100)
	num.mark <- as.integer(tmp$num.mark)
	segmean_ranges <- RangedData(IRanges(start(tmp), end(tmp)),
				     id=sns,
				     seg.mean=seg.mean,
				     num.mark=num.mark,
				     chrom=tmp$chrom)
	colnames(segmean_ranges)[4] <- "chrArm"
	segmean_ranges$chrom <- chromosomeArmToChromosome(segmean_ranges$chrArm)
	segmean_ranges$pedId <- as.factor(who(segmean_ranges$id))
	segmean_ranges <- segmean_ranges[-grep("CIDR", segmean_ranges$pedId), ]
	segmean_ranges$isDeletion <- as.integer(NA)
	return(segmean_ranges)
}

## moved to mybase
##chromosomeArmToChromosome <- function(x){
##	x <- gsub("p", "", x)
##	x <- gsub("q", "", x)
##	return(as.integer(x))
##}

isParent <- function(sample.name, object){
	sample.name <- as.character(sample.name)
	stopifnot("pedId" %in% varLabels(object))
	index <- match(sample.name, sampleNames(object))
	nm <- object$pedId[index]
	nm == "father" | nm == "mother"
}

getDeletions <- function(){
	stopifnot(exists("outdir"))
	segmean_ranges <- checkExists("segmean_ranges", path=outdir, FUN=getSegMeanRanges)
	chr <- chromosomeArmToChromosome(segmean_ranges$chrom)
	## just first 2 chromosomes right now
	for(i in 1:2){#seq(along=chr)){
		index <- which(chr == i)
		object <- segmean_ranges[index, ]
		save(object, file=file.path(outdir, paste(fname, ".rda", sep="")))
	}
	rule.offspring <- function(mad){
		rule <- -2*mad
		rule[rule > log(1/2)] <- log(1/2)
		return(rule)
	}
	rule.parent <- function(mad) -1.5*mad
	object <- callDeletion(object, lrSet, rule.offspring=rule.offspring, rule.parent=rule.parent)
	table(object$isDeletion)
	uid <- unique(object$id)
	del_ranges <- vector("list", length(uid))
	object$chrom <- chromosomeArmToChromosome(object$chrom)
	for(i in seq_along(uid)){
		if(i %% 10 == 0) cat(i, " ")
		subset_ranges <- object[object$id == uid[i], ]
		del_ranges[[i]] <- reduceRD(subset_ranges, by="isDeletion")[[1]]
##		if(i > 1){
##			del_ranges <- c(del_ranges, tmp)
##		} else del_ranges <- tmp
##		rm(tmp)
	}
	del_ranges <- do.call("c", del_ranges)
}
rule.offspring <- function(mad){
	rule <- -2*mad
	rule[rule > log(1/2)] <- log(1/2)
	return(rule)
}
rule.parent <- function(mad) -1.5*mad

sapply.split <- function(x, indices, FUN, ...){
	##sapply(split(x, FACT), FUN, ...)
	## long way more memory efficient?
	##tmp <- split(x, FACT, FUN, ...)
	res <- rep(NA, length(indices))
	for(i in seq_along(indices)){
		if(i %% 1000 == 0) cat(".")
		j <- indices[[i]]
		res[[i]] <- FUN(x[j], ...)
	}
	return(res)
}

getSegMeans <- function(outdir, CHR){
	fnames <- list.files(outdir, pattern=paste("cbs.segs_chr", CHR, "_batch", sep=""), full.name=TRUE)
	if(length(fnames) == 0) stop(paste("There are no segmentation files for chrom", CHR))
	segmeans <- vector("list", length(fnames))
	for(i in seq_along(segmeans)){
		load(fnames[i])
		rd <- RangedData(IRanges(start=cbs.segs$startRow,
					 end=cbs.segs$endRow),
				 pos.start=cbs.segs$loc.start,
				 pos.end=cbs.segs$loc.end,
				 num.mark=cbs.segs$num.mark,
				 seg.mean=cbs.segs$seg.mean,
				 id=cbs.segs$ID,
				 chrom=cbs.segs$chrom)
		segmeans[[i]] <- rd
	}
	segmean_ranges <- do.call("c", segmeans)
	segmean_ranges <- RangedData(IRanges(segmean_ranges$pos.start,
					     segmean_ranges$pos.end),
				     id=segmean_ranges$id,
				     chrom=segmean_ranges$chrom,
				     num.mark=width(segmean_ranges), ## would be number of markers that were not NAs (I think)
				     seg.mean=segmean_ranges$seg.mean,
				     pedId=who(segmean_ranges$id))
	segmean_ranges
}

getRanges <- function(outdir, pattern, CHR, name){
	fnames <- list.files(outdir, pattern=pattern, full.name=TRUE)
	if(missing(name)) stop("must specify object name")
	if(length(fnames) == 0) stop(paste("There are no segmentation files for chrom", CHR))
	segmeans <- vector("list", length(fnames))
	for(i in seq_along(segmeans)){
		load(fnames[i])
		cbs.segs <- get(name)
		rd <- RangedData(IRanges(start=cbs.segs$startRow,
					 end=cbs.segs$endRow),
				 pos.start=cbs.segs$loc.start,
				 pos.end=cbs.segs$loc.end,
				 num.mark=cbs.segs$num.mark,
				 seg.mean=cbs.segs$seg.mean,
				 id=cbs.segs$ID,
				 chrom=cbs.segs$chrom)
		segmeans[[i]] <- rd
	}
	segmean_ranges <- do.call("c", segmeans)
	if(substr(segmean_ranges$id[1], 1, 1) == "X"){
		segmean_ranges$id <- substr(segmean_ranges$id, 2, 9)
	}
	segmean_ranges <- RangedData(IRanges(segmean_ranges$pos.start,
					     segmean_ranges$pos.end),
				     id=segmean_ranges$id,
				     chrom=segmean_ranges$chrom,
				     num.mark=width(segmean_ranges), ## would be number of markers that were not NAs (I think)
				     seg.mean=segmean_ranges$seg.mean,
				     pedId=who(segmean_ranges$id))
	segmean_ranges
}

excludeRanges <- function(segmeans, lrSet){
	##trace(getSegMeanRanges, browser)
	## Drop 168 WGA samples
##	which.wga <- grep("WGA", lrSet$DNA.Source)
	if(length(which.wga) > 0){
		wga.samples <- sampleNames(lrSet)[which.wga]
		segmeans <- segmeans[!segmeans$id %in% wga.samples, ]
	}
	## Drop samples that are not father, mother, offspring
	unknown.ids <- sampleNames(lrSet)[lrSet$pedId == "?"]
	if(length(unknown.ids) > 0)
		segmeans <- segmeans[!segmeans$id %in% unknown.ids, ]
	## Drop other samples in which DNA quality may be affected
	famids.poordna <- lrSet$family[lrSet$MAD > 0.3]
	famids.poordna <- famids.poordna[!is.na(famids.poordna)]
	ids.poordna <- sampleNames(lrSet)[lrSet$family %in% famids.poordna]
	segmeans <- segmeans[!segmeans$id %in% ids.poordna, ]
	segmeans$family <- substr(segmeans$id, 1, 5)
	segmeans$pedId <- who(segmeans$id)
	## remove samples that are not part of a true, or no longer part of a trio as a result of the above filter
	trios <- split(who(segmeans$id), segmeans$family)
	trios <- lapply(trios, unique)
	trios <- trios[sapply(trios, length) == 3]
	family.inTrio <- names(trios)
	segmeans <- segmeans[segmeans$family %in% family.inTrio, ]
	return(segmeans)
}

callDeletion <- function(segmeans, lset, parent.rule, offspring.rule){
	## call for parents
	is.deletion <- rep(NA, nrow(segmeans))
	index <- which(segmeans$pedId == "father" | segmeans$pedId == "mother")
	is.deletion[index] <- ifelse(segmeans$seg.mean[index] < parent.rule(), TRUE, FALSE)
	segmeans$is.deletion <- is.deletion
	if(!"pedId" %in% colnames(segmeans)){
		segmeans$pedId <- who(segmeans$id)
	}
	offspring.index <- which(segmeans$pedId == "offspring")
	offspring.ids <- segmeans$id[offspring.index]
	uids <- unique(offspring.ids)
	offspring.ids <- split(offspring.ids, offspring.ids)
	offspring.ids <- offspring.ids[match(uids, names(offspring.ids))]
	nn <- sapply(offspring.ids, length)

	##uids <- unique(segmeans$id[segmeans$pedId == "offspring"])
	mads <- lset$MAD
	names(mads) <- sampleNames(lset)
	mads <- mads[match(uids, names(mads))]
	stopifnot(identical(names(mads), uids))
	mads <- rep(mads, nn)
	thr <- offspring.rule(mads)
	is.deletion[offspring.index] <- ifelse(segmeans$seg.mean[offspring.index] < thr, TRUE, FALSE)

	## check the BAF if log R ratio is > -1
	index.deletion <- which(is.deletion & segmeans$seg.mean > -1 & segmeans$pedId == "offspring")  ## probably not a homozygous deletion
	fd.ir <- IRanges(position(lset)-12, position(lset)+12)
	deletion.RD <- segmeans[index.deletion, ]
	deletion.ir <- IRanges(start(deletion.RD), end(deletion.RD))
	## for each deletion range, find the indices of lset in the range
	tmp=matchMatrix(findOverlaps(deletion.ir, fd.ir))
	## split by the query index
	## alternatively segment the b allele frequency....
	index.list <- split(tmp[, "subject"], tmp[, "query"])
	sample.index <- match(deletion.RD$id, sampleNames(lset))
	pHet <- rep(NA, length(sample.index))
	for(i in seq_along(index.list)){
		b <- baf(lset)[index.list[[i]], sample.index[i]]
		pHet[i] <- mean(b > 0.2 & b < 0.8, na.rm=TRUE)
	}
	reverseCall <- ifelse(pHet < 0.1, TRUE, FALSE)
	is.deletion[index.deletion] <- reverseCall
	return(is.deletion)
}

##callDenovo <- function(query.list, subject.list, overlapPercentage=0.25){
##	for(i in seq_len(length(query.list))){
##		if(i %% 100 == 0) cat(".")
##		query.ir <- IRanges(start(query.list[[i]]), end(query.list[[i]]))
##		subj.ir <- IRanges(start(subject.list[[i]]), end(subject.list[[i]]))
##		subj.ir <- subj.ir[subject.list[[i]]$is.deletion, ]
##		cnt <- rep(NA, length(query.ir))
##		for(j in 1:length(query.ir)){
##			query <- query.ir[j, ]
##			w <- width(query)
##			if(length(subj.ir) > 0){
##				cnt[j] <- countOverlaps(query, subj.ir, minoverlap=overlapPercentage*w)
##			} else  cnt[j] <- 0
##		}
##		query.list[[i]]$number.overlap <- cnt
##		query.list[[i]]$is.denovo <- cnt==0
##	}
##	return(query.list)
##}

callDenovo <- function(offspring.id, deletion.ranges, epsilon=1e3){
	offspring.index <- which(deletion.ranges$id == offspring.id)
	tmp <- deletion.ranges[offspring.index, ]
	subject.ranges <- IRanges(start(tmp), end(tmp))
	cnt <- countOverlaps(disjoint.ranges, subject.ranges, minoverlap=2L) ## must share start and stop --> minimum of 2
	off.ranges <- disjoint.ranges[cnt > 0, ]

	## do the same thing for the parents
	family.of.offspring <- substr(offspring.id, 1,5)
	parents.of.offspring <- paste(family.of.offspring, c("_02", "_03"), sep="")
	parent.index <- which(substr(deletion.ranges$id, 1, 8) %in% parents.of.offspring)
	tmp <- deletion.ranges[parent.index, ]
	subject.ranges <- IRanges(start(tmp), end(tmp))
	cnt <- countOverlaps(disjoint.ranges, subject.ranges, minoverlap=2L)
	parent.ranges <- disjoint.ranges[cnt > 0, ]

	## extend parent.ranges by epsilon
	parent.ranges <- resize(parent.ranges, width=width(parent.ranges)+2*epsilon, fix="center") ## extends 1kb in each direction

	cnt <- countOverlaps(off.ranges, parent.ranges, minoverlap=2L)
	denovo.ranges <- off.ranges[cnt == 0, ]
	return(denovo.ranges)
}

pBelow <- function(denovo.list, lset){
	for(i in 1:NROW(denovo.list)){
		if(i %% 100 == 0) cat(".")
		query <- denovo.list[[i]]
		father.name <- paste(query$family[1], "_03", sep="")
		mother.name <- paste(query$family[1], "_02", sep="")
		father.index <- grep(father.name, sampleNames(lset))
		mother.index <- grep(mother.name, sampleNames(lset))
		for(j in 1:nrow(query)){
			query2 <- query[j, ]
			p.f <- mean(logR(lset)[start(query2):end(query2), father.index] < parent.rule(), na.rm=TRUE)
			p.m <- mean(logR(lset)[start(query2):end(query2), mother.index] < parent.rule(), na.rm=TRUE)
			query2$propLow.Father <- p.f
			query2$propLow.Mother <- p.m
		}
		denovo.list[[i]] <- query2
	}
	return(denovo.list)
}

## moved to mybase
##myOverlaps <- function(dset.denovo){
##	subj.ir <- IRanges(start(dset.denovo), end(dset.denovo))
##	subj.sns <- dset.denovo$id
##	ix <- order(start(subj.ir))
##	subj.ir <- subj.ir[ix]
##	subj.sns <- subj.sns[ix]
##	brks <- unique(c(start(subj.ir), end(subj.ir)))
##	brks <- brks[order(brks)]
##	ir <- matrix(NA, length(brks)-1, 2)
##	ir[, 1] <- brks[-length(brks)]
##	ir[, 2] <- brks[-1]
##	disj.ir <- IRanges(ir[,1], ir[,2])
##	matching.subj <- vector("list", length(disj.ir))
##	for(i in seq_along(matching.subj)){
##		query.ir <- disj.ir[i, ]
##		w <- width(query.ir)
##		## require overlap to be the length of the interval
##		tmp <- findOverlaps(query.ir, subj.ir, minoverlap=w)
##		matching.subj[[i]] <- matchMatrix(tmp)[, "subject"]
##		##freq[i] <- countOverlaps(query.ir, subj.ir)
##	}
##	disjoint.rd <- RangedData(disj.ir, freq=sapply(matching.subj, length))
##	sns <- lapply(matching.subj, function(i, subj.sns) subj.sns[i], subj.sns)
##	res <- list(disjoint.rd=disjoint.rd, sns=sns)
##	return(res)
##}

constructBeadStudioSet <- function(path, sns, outdir){
}

##plotRange <- function(i, disjoint.rd, binSet, FRAME, samples.altered,
##		      x=c("index", "position"), add.cytoband=TRUE, lset, ...){
##	CHR <- disjoint.rd$chrom[i]
##	rd <- disjoint.rd[i, ]
##	marker.index <- start(rd):end(rd)
##	marker.index <- window(marker.index, FRAME)
##	marker.index <- marker.index[chromosome(binSet)[marker.index] == CHR]
##	if(x[1] == "position"){
##		x <- position(binSet)[marker.index]
##	} else{
##		x <- seq_along(marker.index)
##	}
##	sns <- sampleNames(binSet)
##	if(is(samples.altered, "list")){
##		sample.index <- match(samples.altered[[i]], sns)
##	} else sample.index <- match(samples.altered, sns)
##	y <- loess.residuals(binSet)[marker.index, sample.index]
##	y.raw <- logR(binSet)[marker.index, sample.index]
##	y[y < ylim[1]] <- ylim[1]
##	y[y > ylim[2]] <- ylim[2]
##	y.raw[y.raw < ylim[1]] <- ylim[1]
##	y.raw[y.raw > ylim[2]] <- ylim[2]
##	for(k in seq_along(sample.index)){
##		plot(x, y[, k], pch=21, col="grey60", cex=0.5,xaxt="n", xlab="Mb", ylim=ylim, ...)
##		abline(h=c(THR1, THR2), col="blue", lty=2)
##		tmp <- cbs.ir[cbs.ir$id == colnames(y)[k] & cbs.ir$chrom == CHR, ]
##		sample.segs <- RangedData(IRanges(tmp$loc.start,
##						  tmp$loc.end), seg.mean=tmp$seg.mean)
##		segments(sample.segs, strict=F, lwd=2)
##		legend("topright", legend=paste("grade:", binSet$PanIN.Grade[sample.index[k]]), bty="o",
##		       bg="white")
##		legend("topleft", legend=paste("id:", binSet$individual.id[sample.index[k]]), bty="o",
##		       bg="white")
##		if(add.cytoband){
##			require(SNPchip)
##			data(chromosomeAnnotation)
##			chr.ann <- chromosomeAnnotation[CHR, 1:2]
##			polygon(x=c(chr.ann, rev(chr.ann)),
##				y=c(ylim[1], ylim[1], ylim[2], ylim[2]), col="bisque")
##		}
##		## plot ballele freq.
##		index <- which(chromosome(lset)==CHR & position(lset) >= min(x) & position(lset) <= max(x))
##		ba <- baf(lset)[index, sample.index]
##		plot(position(lset)[index], ba[, k], pch=".", col="grey60")
##	}
##	points(rd$loc.start, ylim[1], pch=24, bg="red", cex=3)
##	points(rd$loc.end, ylim[1], pch=24, bg="red", cex=3)
##	points(rd$loc.start, 0, pch=24, bg="red", cex=3)
##	points(rd$loc.end, 0, pch=24, bg="red", cex=3)
##	at <- pretty(x, n=8)
##	axis(1, at=at, labels=at/1e6, outer=T)
##	mtext(paste("Chr", unique(chromosome(binSet)[marker.index])), 3, outer=TRUE)
##	return()
##}


setMethod("chromosome", "GRanges", function(object) {
	as.integer(sapply(as.character(seqnames(object)), function(x) strsplit(x, "chr")[[1]][2]))
})

featuresInRange <- function(object, range, FRAME=0, FRAME.LEFT, FRAME.RIGHT){
	if(missing(FRAME.LEFT)) FRAME.LEFT <- FRAME
	if(missing(FRAME.RIGHT)) FRAME.RIGHT <- FRAME
 	stopifnot(length(range)==1)
	if(is(range, "GRanges")) CHR <- chromosome(range) else CHR <- range$chrom
	if(FRAME.LEFT > 0 | FRAME.RIGHT > 0){
		require(SNPchip)
		data(chromosomeAnnotation)
		size <- chromosomeAnnotation[CHR, "chromosomeSize"]
		start(range) <- max(start(range) - FRAME.LEFT, 1)
		end(range) <- end(range) + FRAME.RIGHT  ## need to look up chromosome annotation
##		end(range) <- min(end(range), size)
	}
	which(position(object) >= start(range) & position(object) <= end(range) & chromosome(object) == CHR)
}

getFamily <- function(object) substr(sampleNames(object), 1, 5)

triosInRange <- function(object, ## LogRatioSet or something similar
			 range){ ## genomicRanges
	stopifnot(length(range)==1)
	if(is(range, "GRanges")){
		sample.index <- denovoIndicesInRange(range)
	} else {
		sample.index <- match(range$id, sampleNames(object))
	}
	family <- substr(sampleNames(object)[sample.index], 1, 5)
	sampleNames(object)[which(getFamily(object) %in% family)]
}

framePositionIndex <- function(object,  ##LogRatioSet or something similar
			       FRAME){  ##basepairs
		min.pos <- min(pos)-WINDOW.SIZE
	max.pos <- max(pos)+WINDOW.SIZE
}

denovoIndicesInRange <- function(range){
	sample.index <- elementMetadata(range)[, "denovoSamples"]
	as.integer(strsplit(sample.index, ",")[[1]])
}

plotRange <- function(sampleName,    ## names of samples to plot
		      segmentation,  ## the segmentation for the trio
		      lset,          ## LogRatioSet or something similar
		      add.cytoband=TRUE,
		      range,
		      ylim,
		      THR, strict=FALSE, ...){
	stopifnot(length(sampleName) == 1)
	stopifnot(length(range) == 1)
	j <- match(sampleName, sampleNames(lset))
	cn <- copyNumber(lset)[, j]
	if(!missing(ylim)){
		##ylim <- list(...)[["ylim"]]
		cn[cn < ylim[1]] <- ylim[1]
		cn[cn > ylim[2]] <- ylim[2]
		segmentation$seg.mean[segmentation$seg.mean < ylim[1]] <- ylim[1]
		segmentation$seg.mean[segmentation$seg.mean > ylim[2]] <- ylim[2]
	} else ylim <- c(-2,1)
	x <- position(lset)
	plot(x, cn, ylim=ylim, ...)
	abline(h=THR, col="blue", lty=2)
	abline(v=c(start(range), end(range)), lty=2)

	segs <- segmentation[substr(segmentation$id, 1, 8) %in% sampleName, ]
	segments(segs, strict=strict, lwd=2)
	legend("bottomleft", legend=paste("MAD:", round(lset$MAD[j], 3)), bg="white")
	legend("bottomright", legend=who(sampleName),  bty="n", cex=0.8)

	b <- baf(lset)[, j]
	plot(x, b, ylim=c(0,1), ...)
	abline(v=c(start(range), end(range)), lty=2)
	legend("bottomright", legend=who(sampleName),  bty="n", cex=0.8)
}

plotRange2 <- function(sampleName,    ## names of samples to plot
		      segmentation,  ## the segmentation for the trio
		      lset,          ## LogRatioSet or something similar
		      add.cytoband=TRUE,
		      range,
		      ylim,
		      THR, strict=FALSE, ...){
	stopifnot(length(sampleName) == 1)
	stopifnot(length(range) == 1)
	j <- match(sampleName, sampleNames(lset))
	cn <- copyNumber(lset)[, j]
	if(!missing(ylim)){
		##ylim <- list(...)[["ylim"]]
		cn[cn < ylim[1]] <- ylim[1]
		cn[cn > ylim[2]] <- ylim[2]
		segmentation$seg.mean[segmentation$seg.mean < ylim[1]] <- ylim[1]
		segmentation$seg.mean[segmentation$seg.mean > ylim[2]] <- ylim[2]
	} else ylim <- c(-2,1)
	x <- position(lset)
	plot(x, cn, ylim=ylim, ...)
	abline(h=THR, col="blue", lty=2)
	abline(v=c(start(range), end(range)), lty=2)

	segs <- segmentation[segmentation$id %in% sampleName, ]
	segments(segs, strict=strict, lwd=2)
	legend("bottomleft", legend=paste("MAD:", round(lset$MAD[j], 3)), bg="white")
	legend("bottomright", legend=who(sampleName),  bty="n", cex=0.8)
##	b <- baf(lset)[, j]
##	plot(x, b, ylim=c(0,1), ...)
##	abline(v=c(start(range), end(range)), lty=2)
	legend("bottomright", legend=who(sampleName),  bty="n", cex=0.8)
}

##	tmp <- cbs.ir[cbs.ir$id == colnames(y)[k] & cbs.ir$chrom == CHR, ]
##	sample.segs <- RangedData(IRanges(tmp$loc.start,
##					  tmp$loc.end), seg.mean=tmp$seg.mean)
##	segments(sample.segs, strict=F, lwd=2)
##	legend("topright", legend=paste("grade:", binSet$PanIN.Grade[sample.index[k]]), bty="o",
##	       bg="white")
##	legend("topleft", legend=paste("id:", binSet$individual.id[sample.index[k]]), bty="o",
##	       bg="white")
##	if(add.cytoband){
##		require(SNPchip)
##		data(chromosomeAnnotation)
##		chr.ann <- chromosomeAnnotation[CHR, 1:2]
##		polygon(x=c(chr.ann, rev(chr.ann)),
##			y=c(ylim[1], ylim[1], ylim[2], ylim[2]), col="bisque")
##	}
##		## plot ballele freq.
##		index <- which(chromosome(lset)==CHR & position(lset) >= min(x) & position(lset) <= max(x))
##		ba <- baf(lset)[index, sample.index]
##		plot(position(lset)[index], ba[, k], pch=".", col="grey60")
##	}
##	points(rd$loc.start, ylim[1], pch=24, bg="red", cex=3)
##	points(rd$loc.end, ylim[1], pch=24, bg="red", cex=3)
##	points(rd$loc.start, 0, pch=24, bg="red", cex=3)
##	points(rd$loc.end, 0, pch=24, bg="red", cex=3)
##	at <- pretty(x, n=8)
##	axis(1, at=at, labels=at/1e6, outer=T)
##	mtext(paste("Chr", unique(chromosome(binSet)[marker.index])), 3, outer=TRUE)
##	return()
##}
setMethod("$", "GRanges", function(x, name) {
	eval(substitute(elementMetadata(x)$NAME_ARG, list(NAME_ARG=name)))
})

setMethod("colnames", "GRanges", function(x, do.NULL=TRUE, prefix="col") {
	colnames(elementMetadata(x))
})

setReplaceMethod("$", "GRanges", function(x, name, value) {
	elementMetadata(x)[, name] = value
	x
})




freqOfDisjointRange <- function(ranges){
	res <- vector("list", 22)
	for(CHR in 1:22){
		chr.ranges <- ranges[chromosome(ranges) == CHR, ]
		disjoint.ranges <- disjointRanges(chr.ranges)
		subjects <- IRanges(start(chr.ranges),end(chr.ranges))
		cnt <- countOverlaps(disjoint.ranges, subjects, minoverlap=2L)
		disjoint.ranges <- disjoint.ranges[ cnt > 0, ]
		##chr.ranges <- RangedData(chr.ranges, id=ranges$id)
		ranges.ir <- IRanges(start(chr.ranges), end(chr.ranges))

		mm <- matchMatrix(findOverlaps(disjoint.ranges, ranges.ir, minoverlap=2L))
		query.index <- split(1:nrow(mm), mm[, "query"])  ## one disjoint range can overlap many denovo ranges (from different samples)
		median.coverage <- median.size <- rep(NA, length(query.index))      ##  - splitting on the query range groups all of the denovo events that correspond to each disjoint range
		matching.samples <- rep(NA, length(query.index))
		for(j in seq_along(query.index)){
			matching.index <- mm[query.index[[j]], "subject"]
			## list of samples with a deletion in this area.
			tmp <- chr.ranges[matching.index, ]
			matching.samples[j] <- paste(tmp$id, collapse=",")
			segment.sizes <- width(tmp)
			segment.coverage <- tmp$nmarkers##[ii]
			median.size[j] <- median(segment.sizes)
			median.coverage[j] <- median(segment.coverage)
		}
		region <- c(0, cumsum(abs(diff(cnt != 0))))
		chr <- rep(paste("chr", CHR, sep=""), length(disjoint.ranges))
		gr <- GRanges(seqnames=chr,
			      ranges=IRanges(start(disjoint.ranges), end(disjoint.ranges)),
			      freq=cnt[cnt > 0],
			      denovo.samples=matching.samples,
			      median.size=median.size,
			      median.coverage=median.coverage,
			      region=region[cnt>0])
		res[[CHR]] <- gr
	}
	gr <- do.call("c", res)
	return(gr)
}

plotRange4 <- function(ranges, range.index, minDistanceSet, bsSet, outdir){
	mset <- constructTrioSetFromRanges(ranges, minDistanceSet, bsSet)
	CHR <- unique(chromosome(mset))
	ranges.md <- getRanges(outdir,
			       pattern=paste("md.segs.chr", CHR, "_batch", sep=""),
			       name="md.segs", CHR=CHR)
	segmean_ranges <- getSegMeans(outdir, CHR=CHR)
	segmean_ranges$family <- substr(segmean_ranges$id, 1, 5)
	segmean_ranges <- segmean_ranges[segmean_ranges$family %in% sampleNames(mset), ]
	plotSegs(index=range.index,
		 ranges1=ranges,
		 ranges.md=ranges.md,
		 ranges2=segmean_ranges,
		 mset=mset)
}

plotRangeWrapper <- function(i, chrSet, Ranges, segmeans, FRAME, K){
	range.index <- i
	CHR <- chromosome(Ranges)[range.index]
##	if(exists("chrset")){
	if(length(unique(chromosome(chrSet))) > 1){
		load.it <- TRUE
	} else{
		load.it <- ifelse(unique(chromosome(chrSet)) != CHR, TRUE, FALSE)
	}
	tmp <- paste(Ranges$denovo.samples[range.index], collapse=",")
	denovo.samples <- unique(strsplit(tmp, ",")[[1]])
	denovo.families <- substr(denovo.samples, 1, 5)
	if(load.it){
		marker.index <- which(chromosome(chrSet) == CHR)
		all.families <- substr(sampleNames(chrSet), 1, 5)
		j <- which(all.families %in% denovo.families)
		chrSet <- new("LogRatioSet",
			      logRRatio=as.matrix(logR(chrSet)[marker.index, j]),
			      BAF=as.matrix(baf(chrSet)[marker.index, j]),
			      featureData=featureData(chrSet)[marker.index, ],
			      phenoData=phenoData(chrSet)[j, ],
			      annotation=annotation(chrSet))
	}
	if(missing(segmeans)){
		segmeans <- getSegMeans(outdir, CHR=CHR)
	}
	segmeans <- segmeans[substr(segmeans$id, 1, 5) %in% denovo.families, ]
	i <- featuresInRange(chrSet, Ranges[range.index, ], FRAME=FRAME)
	sns <- strsplit(elementMetadata(Ranges)[range.index, "denovo.samples"], ",")[[1]]
##	par(ask=TRUE)
	if(missing(K)) K <- seq_along(sns)
	for(k in K){
		offspring.name <- sns[k]
##		offspring.name <- substr(sns[k], 1, 8)
		offspring.family <- substr(offspring.name, 1, 5)
		parents <- paste(offspring.family, c("_02", "_03"), sep="")
		father.name <- parents[2]
		mother.name <- parents[1]
		jj <- match(c(father.name, mother.name, offspring.name), substr(sampleNames(chrSet), 1, 8))
		lset <- chrSet[i, jj]
		father.name <- sampleNames(lset)[1]; mother.name <- sampleNames(lset)[2]; offspring.name <- sampleNames(lset)[3]
		segmentation <- segmean_ranges[segmean_ranges$id %in% c(offspring.name, mother.name, father.name), ]

		plotRange(father.name, segmentation, lset, range=Ranges[range.index, ],
			  THR=log(1.8/2), ylim=c(-1.5, 0.5),
			  pch=21, col="grey60", cex=0.6, strict=F,
			  xaxt="n")
		plotRange(mother.name, segmentation, lset, range=Ranges[range.index, ],
			  THR=log(1.8/2), ylim=c(-1.5,0.5),
			  pch=21, col="grey60", cex=0.6,
			  xaxt="n")
		THR <- offspring.rule(lset$MAD[match(offspring.name, sampleNames(lset))])
		plotRange(offspring.name, segmentation, lset, range=Ranges[range.index, ],
			  THR=THR, ylim=c(-1.5, 0.5),
			  pch=21, col="grey60", cex=0.6,
			  xaxt="n")
		at <- pretty(range(position(lset)), 8)
		axis(1, at=at, labels=at/1e6)
		mtext(paste("Chr: ", chromosome(Ranges)[range.index], ", offspring ", k , " of ", length(sns), sep=""), 3, outer=T)
	}
	return(list(chrSet=chrSet, segmeans=segmeans))
}

plotRangeWrapper2 <- function(i, chrSet, Ranges, segmeans, FRAME, K,
			      minDistanceSet, distance.ranges){
	sampleNames(chrSet) <- substr(sampleNames(chrSet), 1, 8)
	range.index <- i[[1]]
	rm(i)
	CHR <- chromosome(Ranges)[range.index]
##	if(exists("chrset")){
	if(length(unique(chromosome(chrSet))) > 1){
		stop("length(unique(chromosome)) > 1")
		load.it <- TRUE
	} else{
		load.it <- ifelse(unique(chromosome(chrSet)) != CHR, TRUE, FALSE)
	}
	tmp <- paste(Ranges$denovo.samples[range.index], collapse=",")
	denovo.samples <- unique(strsplit(tmp, ",")[[1]])
	denovo.families <- substr(denovo.samples, 1, 5)
##	if(load.it){
##		marker.index <- which(chromosome(chrSet) == CHR)
##		all.families <- substr(sampleNames(chrSet), 1, 5)
##		j <- which(all.families %in% denovo.families)
##		chrSet <- new("LogRatioSet",
##			      logRRatio=as.matrix(logR(chrSet)[marker.index, j]),
##			      BAF=as.matrix(baf(chrSet)[marker.index, j]),
##			      featureData=featureData(chrSet)[marker.index, ],
##			      phenoData=phenoData(chrSet)[j, ],
##			      annotation=annotation(chrSet))
##	}
	if(missing(segmeans)){
		segmeans <- getSegMeans(outdir, CHR=CHR)
	}
	segmeans <- segmeans[substr(segmeans$id, 1, 5) %in% denovo.families, ]
	i <- featuresInRange(chrSet, Ranges[range.index, ], FRAME=FRAME)
	sns <- strsplit(elementMetadata(Ranges)[range.index, "denovo.samples"], ",")[[1]]
##	par(ask=TRUE)
	if(missing(K)) K <- seq_along(sns)
	for(k in K){
		offspring.name <- substr(sns[k], 1, 8)
##		offspring.name <- substr(sns[k], 1, 8)
		offspring.family <- substr(offspring.name, 1, 5)
		parents <- paste(offspring.family, c("_02", "_03"), sep="")
		father.name <- parents[2]
		mother.name <- parents[1]
		jj <- match(c(father.name, mother.name, offspring.name), substr(sampleNames(chrSet), 1, 8))
		lset <- chrSet[i, jj]
		iii <- match(featureNames(lset), featureNames(minDistanceSet))
		jjj <- match(offspring.name, sampleNames(minDistanceSet))
		mset <- minDistanceSet[iii, jjj]##copyNumber(minDistanceSet)[iii, jjj]
		father.name <- sampleNames(lset)[1]; mother.name <- sampleNames(lset)[2]; offspring.name <- sampleNames(lset)[3]
		segmentation <- segmean_ranges[substr(segmean_ranges$id, 1, 8) %in% c(offspring.name, mother.name, father.name), ]
		plotRange(father.name, segmentation, lset, range=Ranges[range.index, ],
			  THR=log(1.8/2), ylim=c(-1.5, 0.5),
			  pch=21, col="grey60", cex=0.6, strict=F,
			  xaxt="n")
		plotRange(mother.name, segmentation, lset, range=Ranges[range.index, ],
			  THR=log(1.8/2), ylim=c(-1.5,0.5),
			  pch=21, col="grey60", cex=0.6,
			  xaxt="n")
		offspring.rule2 <- function(MAD) ifelse(-1.5*MAD < log(1.5/2), -1.5*MAD, log(1.5/2))
		THR <- offspring.rule2(lset$MAD[match(offspring.name, sampleNames(lset))])
		plotRange(offspring.name, segmentation, lset, range=Ranges[range.index, ],
			  THR=THR, ylim=c(-1.5, 0.5),
			  pch=21, col="grey60", cex=0.6,
			  xaxt="n")
		plotRange2(offspring.name, distance.ranges, mset, range=Ranges[range.index, ],
			   THR=THR, ylim=c(-0.5, 1.5),
			   pch=21, col="grey60", cex=0.6,
			   xaxt="n")
		at <- pretty(range(position(lset)), 8)
		axis(1, at=at, labels=at/1e6)
		mtext(paste("Chr: ", chromosome(Ranges)[range.index], ", offspring ", k , " of ", length(sns), sep=""), 3, outer=T)
	}
	##return(list(chrSet=chrSet, segmeans=segmeans))
	NULL
}

plotSegs <- function(index, ranges1, ranges.md, ranges2, mset, FRAME, strict=FALSE, ylim=c(-1.5, 0.5), ...){
	require(SNPchip)
	data(chromosomeAnnotation)
	ranges2$id <- substr(ranges2$id, 1, 8)
	ranges.md <- ranges.md[ranges.md$id == ranges1$id[index], ]
	range.index <- index[[1]]
	ranges1$family <- substr(ranges1$id, 1, 5)
	CHR <- ranges1$chrom[index]
	chrAnn <- chromosomeAnnotation[CHR, ]
	this.range <- ranges1[index, ]
	ranges1 <- ranges1[ranges1$id %in% this.range$id, ]
	ranges2 <- ranges2[ranges2$family %in% this.range$family, ]
	ranges2$seg.mean[ranges2$seg.mean < ylim[1]] <- ylim[1]
	ranges2$seg.mean[ranges2$seg.mean > ylim[2]] <- ylim[2]
	## goal is to have the feature cover approx 5% of the plot
	if(missing(FRAME)){
		w <- width(this.range)
		FRAME <- w/0.05  * 1/2
	}
	i <- featuresInRange(mset, this.range, FRAME=FRAME)
	family.name <- this.range$id
	j <- match(this.range$family, sampleNames(mset))
	lset <- mset[i, j]
	ranges.F <- ranges2[grep("_03", ranges2$id), ]
	ranges.M <- ranges2[grep("_02", ranges2$id), ]
	ranges.O <- ranges2[grep("_01", ranges2$id), ]
	xx <- c(chrAnn[1:2], chrAnn[2:1])
	yy <- c(-1.2, -1.2, 0.4, 0.4)
	yyc <- c(0.1, 0.1, 0.9, 0.9)
	x <- position(lset)
	y <- logR.F(lset)
	y[y < ylim[1]] <- ylim[1]
	y[y > ylim[2]] <- ylim[2]
	##FATHER
	plot(x, y, pch=21, col="blue", cex=0.6, xaxt="n", ylim=ylim)
	abline(v=c(start(this.range), end(this.range)), lty=2)
	segments(ranges.F, strict=strict)
	polygon(xx, yy, col="bisque")
	plot(x, baf.F(lset), pch=21, col="red", cex=0.6, xaxt="n", ylim=c(0,1))
	abline(v=c(start(this.range), end(this.range)), lty=2)
	polygon(xx, yyc, col="bisque")
	##MOTHER
	y <- logR.M(lset)
	y[y < ylim[1]] <- ylim[1]
	y[y > ylim[2]] <- ylim[2]
	plot(x, y, pch=21, col="blue", cex=0.6, xaxt="n", ylim=ylim)
	abline(v=c(start(this.range), end(this.range)), lty=2)
	segments(ranges.M, strict=strict)
	polygon(xx, yy, col="bisque")
	plot(x, baf.M(lset), pch=21, col="red", cex=0.6, xaxt="n", ylim=c(0,1))
	abline(v=c(start(this.range), end(this.range)), lty=2)
	polygon(xx, yyc, col="bisque")
	##OFFSPRING
	y <- logR.O(lset)
	y[y < ylim[1]] <- ylim[1]
	y[y > ylim[2]] <- ylim[2]
	plot(x, y, pch=21, col="blue", cex=0.6, xaxt="n", ylim=ylim)
	abline(v=c(start(this.range), end(this.range)), lty=2)
	segments(ranges.O, strict=strict)
	legend("topleft", legend=paste("MAD:", round(lset$MAD,3)), bty="n")
	polygon(xx, yy, col="bisque")
	plot(x, baf.O(lset), pch=21, col="red", cex=0.6, xaxt="n", ylim=c(0,1))
	abline(v=c(start(this.range), end(this.range)), lty=2)
	polygon(xx, yyc, col="bisque")
	yym <- c(-0.4, -0.4, 0.9, 0.9)
	y <- mindist(lset)
	y[y < -0.5] <- -0.5
	y[y > 1] <- 1
	plot(x, y, pch=21, col="grey50", cex=0.6, xaxt="n", ylim=c(-0.5,1))
	abline(v=c(start(this.range), end(this.range)), lty=2)
	segments(ranges.md, strict=strict)
	polygon(xx, yym, col="bisque")
	at <- pretty(range(position(lset)), 8)
	axis(1, at=at, labels=at/1e6)
	mtext(paste("Family ", substr(this.range$id, 1,5), "; Chr: ", CHR, sep=""), 3, outer=T)
}

plotSegs2 <- function(index, ranges1, ranges.md, ranges2,
		      mset, FRAME, strict=FALSE, ylim=c(-1.5, 0.5), ...){
	require(SNPchip)
	data(chromosomeAnnotation)
	ranges2$id <- substr(ranges2$id, 1, 8)
	ranges.md <- ranges.md[ranges.md$id == paste(sampleNames(mset), "_01", sep=""), ]
	##ranges.md <- ranges.md[ranges.md$id == ranges1$id[index], ]
	range.index <- index[[1]]
	ranges1$family <- substr(ranges1$id, 1, 5)
	CHR <- ranges1$chrom[index]
	chrAnn <- chromosomeAnnotation[CHR, ]
	this.range <- ranges1[index, ]
	##ranges1 <- ranges1[ranges1$id %in% this.range$id, ]
	##ranges2 <- ranges2[ranges2$family %in% substr(sampleNames(mset), 1, 5), ]
	ranges2$seg.mean[ranges2$seg.mean < ylim[1]] <- ylim[1]
	ranges2$seg.mean[ranges2$seg.mean > ylim[2]] <- ylim[2]
	## goal is to have the feature cover approx 5% of the plot
	if(missing(FRAME)){
		w <- width(this.range)
		FRAME <- w/0.05  * 1/2
	}
	i <- featuresInRange(mset, this.range, FRAME=FRAME)
	##family.name <- this.range$id
	family.name <- substr(sampleNames(mset), 1, 5)
	##j <- match(family.name, sampleNames(mset))
	lset <- mset[i, ]
	ranges.F <- ranges2[grep("_03", ranges2$id), ]
	ranges.M <- ranges2[grep("_02", ranges2$id), ]
	ranges.O <- ranges2[grep("_01", ranges2$id), ]
	xx <- c(chrAnn[1:2], chrAnn[2:1])
	yy <- c(-1.2, -1.2, 0.4, 0.4)
	yyc <- c(0.1, 0.1, 0.9, 0.9)
	x <- position(lset)
	y <- logR.F(lset)
	y[y < ylim[1]] <- ylim[1]
	y[y > ylim[2]] <- ylim[2]
	##FATHER
	plot(x, y, pch=21, col="blue", cex=0.6, xaxt="n", ylim=ylim)
	abline(v=c(start(this.range), end(this.range)), lty=2)
	segments(ranges.F, strict=strict)
	polygon(xx, yy, col="bisque")
	plot(x, baf.F(lset), pch=21, col="red", cex=0.6, xaxt="n", ylim=c(0,1))
	abline(v=c(start(this.range), end(this.range)), lty=2)
	polygon(xx, yyc, col="bisque")
	##MOTHER
	y <- logR.M(lset)
	y[y < ylim[1]] <- ylim[1]
	y[y > ylim[2]] <- ylim[2]
	plot(x, y, pch=21, col="blue", cex=0.6, xaxt="n", ylim=ylim)
	abline(v=c(start(this.range), end(this.range)), lty=2)
	segments(ranges.M, strict=strict)
	polygon(xx, yy, col="bisque")
	plot(x, baf.M(lset), pch=21, col="red", cex=0.6, xaxt="n", ylim=c(0,1))
	abline(v=c(start(this.range), end(this.range)), lty=2)
	polygon(xx, yyc, col="bisque")
	##OFFSPRING
	y <- logR.O(lset)
	y[y < ylim[1]] <- ylim[1]
	y[y > ylim[2]] <- ylim[2]
	plot(x, y, pch=21, col="blue", cex=0.6, xaxt="n", ylim=ylim)
	abline(v=c(start(this.range), end(this.range)), lty=2)
	segments(ranges.O, strict=strict)
	polygon(xx, yy, col="bisque")
	plot(x, baf.O(lset), pch=21, col="red", cex=0.6, xaxt="n", ylim=c(0,1))
	abline(v=c(start(this.range), end(this.range)), lty=2)
	polygon(xx, yyc, col="bisque")
	legend("topleft", legend=paste("MAD:", round(lset$MAD,3)), bty="n")
	yym <- c(-0.4, -0.4, 0.9, 0.9)
	y <- mindist(lset)
	y[y < -0.5] <- -0.5
	y[y > 1] <- 1
	plot(x, y, pch=21, col="grey50", cex=0.6, xaxt="n", ylim=c(-0.5,1))
	abline(v=c(start(this.range), end(this.range)), lty=2)
	segments(ranges.md, strict=strict)
	polygon(xx, yym, col="bisque")
	at <- pretty(range(position(lset)), 8)
	axis(1, at=at, labels=at/1e6)
	mtext(paste("Family ", substr(this.range$id, 1,5), "; Chr: ", CHR, sep=""), 3, outer=T)
}

plotRange3 <- function(sampleName,    ## names of samples to plot
		      segmentation,  ## the segmentation for the trio
		      lset,          ## LogRatioSet or something similar
		      add.cytoband=TRUE,
		      range,
		      ylim,
		      THR, strict=FALSE, ...){
	stopifnot(length(sampleName) == 1)
	stopifnot(length(range) == 1)
	j <- match(sampleName, sampleNames(lset))
	cn <- copyNumber(lset)[, j]
	if(!missing(ylim)){
		##ylim <- list(...)[["ylim"]]
		cn[cn < ylim[1]] <- ylim[1]
		cn[cn > ylim[2]] <- ylim[2]
		segmentation$seg.mean[segmentation$seg.mean < ylim[1]] <- ylim[1]
		segmentation$seg.mean[segmentation$seg.mean > ylim[2]] <- ylim[2]
	} else ylim <- c(-2,1)
	x <- position(lset)
	plot(x, cn, ylim=ylim, ...)
	abline(h=THR, col="blue", lty=2)
	abline(v=c(start(range), end(range)), lty=2)

	segs <- segmentation[substr(segmentation$id, 1, 8) %in% sampleName, ]
	segments(segs, strict=strict, lwd=2)
	legend("bottomleft", legend=paste("MAD:", round(lset$MAD[j], 3)), bg="white")
	legend("bottomright", legend=who(sampleName),  bty="n", cex=0.8)

	b <- baf(lset)[, j]
	plot(x, b, ylim=c(0,1), ...)
	abline(v=c(start(range), end(range)), lty=2)
	legend("bottomright", legend=who(sampleName),  bty="n", cex=0.8)
}


overlapsCentromere <- function(myranges, centromere.ranges, CHR){
	myranges.bak <- myranges
	centromere.ir <- IRanges(start(centromere.ranges)[CHR],
				 end(centromere.ranges)[CHR])
	if(!is(class(myranges), "IRanges")){
		myranges <- IRanges(start(myranges), end(myranges))
	}
	cnt <- countOverlaps(myranges, centromere.ir)
##	if(length(index) > 0) {
##		myranges <- myranges.bak[index, ]
##	} else myranges <- myranges.bak
	return(cnt > 0)
}

myunlist <- function(rdList){
	id <- rep(names(rdList), sapply(rdList, length))
	starts <- unlist(lapply(rdList, start))
	ends <- unlist(lapply(rdList, end))
	denovo.ir <- IRanges(starts, ends)
	denovo.ranges <- RangedData(denovo.ir, id=id)
}


statisticsForRanking <- function(deletion.ranges, disjoint.ranges, CHR){
	deletion.ir <- IRanges(start(deletion.ranges), end(deletion.ranges))
	cnt <- countOverlaps(disjoint.ranges, deletion.ir, minoverlap=2L)
	##disjoint.ranges <- disjoint.ranges[cnt > 0, ]  ## makes it run faster
	mm <- matchMatrix(findOverlaps(disjoint.ranges, deletion.ir, minoverlap=2L))
	query.index <- split(1:nrow(mm), mm[, "query"])  ## one disjoint range can overlap many denovo ranges (from different samples)
	median.proportion.overlap <- median.coverage <- median.size <- rep(NA, length(query.index))      ##  - splitting on the query range groups all of the denovo events that correspond to each disjoint range
	denovo.samples <- rep(NA, length(query.index))
	for(j in seq_along(query.index)){
		## if(j %% 10 == 0) cat(j, " ")
		matching.index <- mm[query.index[[j]], "subject"]
		## list of samples with a deletion in this area.
		tmp <- deletion.ranges[matching.index, ]
		denovo.samples[j] <- paste(substr(tmp$id, 1, 8),  collapse=",")
		##query <- IRanges(unique(start(tmp)), unique(end(tmp)))
		query <- IRanges(start(tmp), end(tmp))
		tmp2 <- deletion.ranges[deletion.ranges$id %in% tmp$id, ]
		subject <- IRanges(start(tmp2), end(tmp2))
		ii <- matchMatrix(findOverlaps(query, subject, minoverlap=2L))[, "subject"]
		segment.sizes <- width(tmp2)[ii]
		segment.coverage <- tmp2$num.mark[ii]
		median.size[j] <- median(segment.sizes)
		median.coverage[j] <- median(segment.coverage)
	}
	region <- c(0, cumsum(abs(diff(cnt != 0))))
	region <- region[cnt > 0]
	chr <- paste("chr", CHR, sep="")
	gr <- GRanges(seqnames=Rle(chr, length(disjoint.ranges)),
		      ranges=IRanges(start(disjoint.ranges),
		      end(disjoint.ranges)),
		      freq=cnt[cnt > 0],
		      denovo.samples=denovo.samples,
		      median.size=median.size,
		      median.coverage=median.coverage,
		      region=region)
}

callDenovoAllTrios <- function(deletion.ranges, disjoint.ranges, epsilon=2){
	offspring.samples <- unique(deletion.ranges$id[deletion.ranges$pedId == "offspring"])
	rdList <- vector("list", length(offspring.samples))
	for(i in seq_along(offspring.samples)){
		##		trace(callDenovo, browser)
		rdList[[i]] <- callDenovo(offspring.samples[i], deletion.ranges, epsilon=epsilon)
	}
	names(rdList) <- offspring.samples
	## denovo.ranges are a subset of the set of disjoint ranges
	denovo.ranges <- myunlist(rdList)
	denovo.ranges
}

callDeletion2 <- function(ranges, trioSet, offspring.rule){
	mads <- trioSet$mad
	offspring.ids <- split(1:nrow(ranges), ranges$id)
	nn <- sapply(offspring.ids, length)
	thr <- rep(offspring.rule(mads), nn)
	is.deletion <- ifelse(ranges$seg.mean > thr, TRUE, FALSE)

	index.deletion <- which(is.deletion & ranges$seg.mean > 0)
	##candidate deletion ranges
	deletion.RD <- ranges[index.deletion, ]
	deletion.ir <- IRanges(start(deletion.RD), end(deletion.RD))

	##feature data ranges
	fD <- fData(trioSet)
	fd.ir <- IRanges(fD$position-12, fD$position+12)
	tmp <- matchMatrix(findOverlaps(deletion.ir, fd.ir))
	B <- assayData(trioSet)[["baf.O"]]
	## split by the query index
	## alternatively segment the b allele frequency....
	index.list <- split(tmp[, "subject"], tmp[, "query"])
	sample.index <- match(substr(deletion.RD$id, 1, 5), sampleNames(trioSet))
	pHet <- rep(NA, length(sample.index))
	for(i in seq_along(index.list)){
		##	b <- baf(lset)[index.list[[i]], sample.index[i]]
		b <- B[index.list[[i]], sample.index[i]]
		pHet[i] <- mean(b > 0.2 & b < 0.8, na.rm=TRUE)
	}
	is.deletion[index.deletion] <- ifelse(pHet < 0.1, TRUE, FALSE)
	return(is.deletion)
}


minDistanceDeletion <- function(ranges, minDistanceSet, offspring.rule, CHR){
	mads <- minDistanceSet$Mad
	offspring.ids <- split(1:nrow(ranges), ranges$id)
	nn <- sapply(offspring.ids, length)
	thr <- rep(offspring.rule(mads), nn)
	is.deletion <- ifelse(ranges$seg.mean > thr, TRUE, FALSE)

	## check the BAF if log R ratio is > -1
##	sampleNames(bsSet) <- substr(sampleNames(bsSet), 1, 8)
	##I <- match(featureNames(minDistanceSet), featureNames(bsSet))

##	I <- which(chromosome(minDistanceSet) == CHR)
##	J <- match(sampleNames(minDistanceSet), sampleNames(bsSet))
##	B <- as.matrix(baf(bsSet)[I, J])
##	fD <- fData(bsSet)[I, ]
##	invisible(open(baf(minDistanceSet)))
##	B <- as.matrix(baf(minDistanceSet)[I, J])
	B <- baf(minDistanceSet)
##	invisible(close(baf(minDistanceSet)))
##	fD <- fData(minDistanceSet)[I, ]
	fD <- fData(minDistanceSet)

	index.deletion <- which(is.deletion & ranges$seg.mean > -1)  ## probably not a homozygous deletion
	fd.ir <- IRanges(fD$position-12, fD$position+12)
	deletion.RD <- ranges[index.deletion, ]
	deletion.ir <- IRanges(start(deletion.RD), end(deletion.RD))
	tmp <- matchMatrix(findOverlaps(deletion.ir, fd.ir))
	## split by the query index
	## alternatively segment the b allele frequency....
	index.list <- split(tmp[, "subject"], tmp[, "query"])
	sample.index <- match(deletion.RD$id, sampleNames(minDistanceSet))
	pHet <- rep(NA, length(sample.index))
	for(i in seq_along(index.list)){
		##	b <- baf(lset)[index.list[[i]], sample.index[i]]
		b <- B[index.list[[i]], sample.index[i]]
		pHet[i] <- mean(b > 0.2 & b < 0.8, na.rm=TRUE)
	}
	is.deletion[index.deletion] <- ifelse(pHet < 0.1, TRUE, FALSE)
	return(is.deletion)
}


readPennCnv <- function(penndir){
	fnames <- list.files(penndir, full.names=TRUE)
	tmp <- read.delim(fnames[1], nrows=5, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	colClasses <- c("integer", "integer", "integer", "integer", "character", "character",
			"character", "character", "character",
			"character")
	rdList <- vector("list", length(fnames))
	for(i in seq_along(fnames)){
		if(i %% 100 == 0) cat('.')
		penn.joint <- read.delim(fnames[i], header=TRUE, sep="\t", stringsAsFactors=FALSE,
					 colClasses=colClasses)
		penn.joint$LengthCNV <- as.integer(gsub(",", "", penn.joint$LengthCNV))
		chr <- paste("chr", penn.joint$Chromosome, sep="")
		gr <- GRanges(seqnames=chr,
			      ranges=IRanges(penn.joint$StartPosition,
			      penn.joint$EndPosition),
			      nmarkers=penn.joint$NumberSNPs,
			      id=penn.joint$ID,
			      triostate=penn.joint$TrioState,
			      pedId=penn.joint$FamilyMember)
		rdList[[i]] <- gr
	}
##	pennJoint <- do.call("c", rdList)
##	return(pennJoint)
	rdList
}

constructTrioSet <- function(minDistanceSet, bsSet, CHR){
	J <- match(sampleNames(minDistanceSet), sampleNames(bsSet))
	I <- which(chromosome(minDistanceSet) == CHR)
	stopifnot(identical(featureNames(minDistanceSet)[I], featureNames(bsSet)[I]))
	sample.names <- substr(sampleNames(minDistanceSet), 1, 5)
	father.names <- paste(sample.names, "03", sep="_")
	mother.names <- paste(sample.names, "02", sep="_")
	father.index <- match(father.names, sampleNames(bsSet))
	mother.index <- match(mother.names, sampleNames(bsSet))
	offspr.index <- match(sampleNames(minDistanceSet),sampleNames(bsSet))
	logR.F <- as.matrix(logR(bsSet)[I, father.index])
	logR.M <- as.matrix(logR(bsSet)[I, mother.index])
	logR.O <- as.matrix(logR(bsSet)[I, offspr.index])
	baf.F <- as.matrix(baf(bsSet)[I, father.index])
	baf.M <- as.matrix(baf(bsSet)[I, mother.index])
	baf.O <- as.matrix(baf(bsSet)[I, offspr.index])
	colnames(logR.F) <- colnames(logR.M) <- colnames(logR.O) <- sample.names
	colnames(baf.F) <- colnames(baf.M) <- colnames(baf.O) <- sample.names
	phenoD <- phenoData(bsSet)[offspr.index, ]
	sampleNames(phenoD) <- sample.names
	mindist <- as.matrix(copyNumber(minDistanceSet)[I, ])
	colnames(mindist) <- sample.names
	mads <- apply(mindist, 2, mad, na.rm=TRUE)
	mset <- new("MultiSet",
		    mindist=mindist,
		    logR.F=logR.F,
		    logR.M=logR.M,
		    logR.O=logR.O,
		    baf.F=baf.F,
		    baf.M=baf.M,
		    baf.O=baf.O,
		    phenoData=phenoD,
		    featureData=featureData(bsSet)[I, ])
	mset$mad <- mads
	rm(logR.F, logR.M, logR.O, baf.F, baf.M, baf.O, mindist); gc()
	mset$mad <- mads
	return(mset)
}


constructTrioSetFromRanges <- function(ranges1, ## top hit ranges
				       ##ranges2, ## entire segmentation
				       minDistanceSet,
				       bsSet,
				       FRAME,
				       xlim,
				       id){
	require(SNPchip)
	data(chromosomeAnnotation)
	## marker indices
	marker.index <- list()
	for(i in 1:nrow(ranges1)){
		CHR <- ranges1$chrom[i]
		chrAnn <- chromosomeAnnotation[CHR, ]
		this.range <- ranges1[i, ]
##		ranges1 <- ranges1[ranges1$id %in% this.range$id, ]
##		ranges2 <- ranges2[ranges2$family %in% this.range$family, ]
##		ranges2$seg.mean[ranges2$seg.mean < ylim[1]] <- ylim[1]
##		ranges2$seg.mean[ranges2$seg.mean > ylim[2]] <- ylim[2]
		## goal is to have the feature cover approx 5% of the plot
		if(missing(xlim)){
			if(missing(FRAME)){
				w <- width(this.range)
				FRAME <- w/0.05  * 1/2
			}
			marker.index[[i]] <- featuresInRange(minDistanceSet, this.range, FRAME=FRAME)
		} else {
			marker.index[[i]] <- which(position(minDistanceSet) >= min(xlim) & position(minDistanceSet) <= max(xlim) & chromosome(minDistanceSet)==CHR)
		}
	}
##	ls <- sapply(marker.index, length)
##	regions <- rep(1:length(ls), ls)
	marker.index <- unique(unlist(marker.index))
	marker.index <- marker.index[order(marker.index)]
	## sample indices
	if(missing(id)){
		sample.index <- list()
		for(i in 1:nrow(ranges1)){
			sample.index[[i]] <- grep(substr(ranges1$id[i], 1, 5), sampleNames(bsSet))
		}
		sample.index <- unique(unlist(sample.index))
	} else {
		sample.index <- match(id, sampleNames(bsSet))
	}
	J <- sample.index
	I <- marker.index

	##J <- match(sampleNames(minDistanceSet), sampleNames(bsSet))
	##I <- which(chromosome(minDistanceSet) == CHR)
	stopifnot(identical(featureNames(minDistanceSet)[I], featureNames(bsSet)[I]))
	##sample.names <- substr(sampleNames(minDistanceSet), 1, 5)
	if(missing(id)){
		sample.names <- unique(substr(ranges1$id, 1, 5))
	} else sample.names <- substr(id, 1, 5)
	father.names <- paste(sample.names, "03", sep="_")
	mother.names <- paste(sample.names, "02", sep="_")
	father.index <- match(father.names, sampleNames(bsSet))
	mother.index <- match(mother.names, sampleNames(bsSet))
	offspr.index <- match(paste(sample.names, "01", sep="_"), sampleNames(bsSet))
	##offspr.index <- match(unique(ranges1$id), sampleNames(bsSet))
	##offspr.index <- match(sampleNames(minDistanceSet),sampleNames(bsSet))

	logR.F <- as.matrix(logR(bsSet)[I, father.index])
	logR.M <- as.matrix(logR(bsSet)[I, mother.index])
	logR.O <- as.matrix(logR(bsSet)[I, offspr.index])
	baf.F <- as.matrix(baf(bsSet)[I, father.index])
	baf.M <- as.matrix(baf(bsSet)[I, mother.index])
	baf.O <- as.matrix(baf(bsSet)[I, offspr.index])
	colnames(logR.F) <- colnames(logR.M) <- colnames(logR.O) <- sample.names
	colnames(baf.F) <- colnames(baf.M) <- colnames(baf.O) <- sample.names
	phenoD <- phenoData(bsSet)[offspr.index, ]
	sampleNames(phenoD) <- sample.names

	offspr.index <- match(paste(sample.names, "01", sep="_"), sampleNames(minDistanceSet))
	mindist <- as.matrix(copyNumber(minDistanceSet)[I, offspr.index])
	colnames(mindist) <- sample.names
	##mads <- apply(mindist, 2, mad, na.rm=TRUE)
	mset <- new("MultiSet",
		    mindist=mindist,
		    logR.F=logR.F,
		    logR.M=logR.M,
		    logR.O=logR.O,
		    baf.F=baf.F,
		    baf.M=baf.M,
		    baf.O=baf.O,
		    phenoData=phenoD,
		    featureData=featureData(bsSet)[I, ])
	index <- match(sampleNames(mset), sampleNames(bsSet))
	mset$mad <- bsSet$MAD[index]
	##mset$mad <- mads
	rm(logR.F, logR.M, logR.O, baf.F, baf.M, baf.O, mindist); gc()
	##mset$mad <- mads
	return(mset)
}

collectAllRangesOfSize <- function(SIZE, bsSet, minDistanceSet,
				   ##offspring.rule,
				   outdir, MIN=1, MAX=4, lambda=0.1){
##	bsSet <- checkExists("bsSet", .path=outdir, .FUN=load)
##	sampleNames(bsSet) <- substr(sampleNames(bsSet), 1, 8)
##	open(baf(bsSet))
##	open(logR(bsSet))
	library(SNPchip)
	data(chromosomeAnnotation)
	centromere.ranges <- GRanges(seqnames=Rle(paste("chr", 1:22, sep=""), rep(1,22)),
				     ranges=IRanges(chromosomeAnnotation[1:22, "centromereStart"],
				     chromosomeAnnotation[1:22, "centromereEnd"]))
##	minDistanceSet <- checkExists("minDistanceSet", .path=outdir, .FUN=load)
##	invisible(open(copyNumber(minDistanceSet)))
	stopifnot(identical(featureNames(bsSet), featureNames(minDistanceSet)))
	deletion.ranges <- vector("list", 22)
	for(CHR in 1:22){
		cat(CHR, "\n")
		ranges <- getRanges(outdir, pattern=paste("md.segs.chr", CHR, "_batch", sep=""),
				    name="md.segs", CHR=CHR)
		ranges2 <- ranges[ranges$seg.mean > 0 & ranges$num.mark >= SIZE, ]
		index <- match(ranges2$id, sampleNames(minDistanceSet))
		mads <- minDistanceSet$Mad[index]

		x <- ranges2$num.mark
##		f <- function(x, lambda, MIN, MAX){
		p <- lambda*exp(-lambda*x)
		##rescale p to have range [1, 4]
		MIN <- 1; MAX <- 4
		b <- 1/(MAX - MIN)
		a <- MIN * b
		numberMads <- ((p-min(p))/(max(p)-min(p)) + a)/b
##		numberMads
##		}
##		thr <- offspring.rule(mads)
##		thr[thr < 0.2] <- 0.2
##		thr <- p.rescaled * mads
		thr <- numberMads * mads
		thr[thr < 0.2] <- 0.2
		ranges2$is.deletion <- ifelse(ranges2$seg.mean >= thr, TRUE, FALSE)
		deletion.RD <- ranges2[ranges2$is.deletion, ]
		deletion.ir <- IRanges(start(deletion.RD), end(deletion.RD))
		fD <- fData(minDistanceSet)[chromosome(minDistanceSet) == CHR, ]
		fd.ir <- IRanges(fD$position-12, fD$position+12)
		tmp <- matchMatrix(findOverlaps(deletion.ir, fd.ir))
		sample.index <- match(deletion.RD$id, sampleNames(bsSet))
		index.list <- split(tmp[, "subject"], tmp[, "query"])
		fns.list <- lapply(index.list, function(i, fns) fns[i], fns=rownames(fD))
		chrom.index <- which(chromosome(bsSet) == CHR)
		B <- as.matrix(baf(bsSet)[chrom.index[unique(as.integer(unlist(index.list)))], sample.index])
		pHet <- rep(NA, length(sample.index))
		for(i in seq_along(index.list)){
			ii <- match(fns.list[[i]], rownames(B))
			b <- B[ii, i]
			pHet[i] <- mean(b > 0.2 & b < 0.8, na.rm=TRUE)
		}
		deletion.RD$pHet <- pHet
		deletion.RD$noCentromere.overlap <- !(overlapsCentromere(deletion.RD, centromere.ranges, CHR))
		deletion.RD <- deletion.RD[order(start(deletion.RD)), ]
		##ranges <- ranges[order(start(ranges), end(ranges)), ]
		deletion.RD$is.deletion <- deletion.RD$pHet < 0.05 & deletion.RD$is.deletion & deletion.RD$noCentromere.overlap
		deletion.ranges[[CHR]] <- deletion.RD
		rm(ranges, ranges2, deletion.RD, B, deletion.ir);gc()
	}
	tmp <- do.call("c", deletion.ranges)
	deletion.ranges <- RangedData(IRanges(start(tmp), end(tmp)),
				      id=tmp$id,
				      chrom=tmp$chrom,
				      num.mark=tmp$num.mark,
				      seg.mean=tmp$seg.mean,
				      pHet=tmp$pHet,
				      is.deletion=tmp$is.deletion,
				      noCentromere.overlap=tmp$noCentromere.overlap)
	return(deletion.ranges)
}

getRefGene <- function(filename="~/Data/Downloads/hg18_refGene.txt"){
	##tmp <- read.delim("~/Downloads/hg18_refGene.txt", nrows=5, header=FALSE)
	##colnames(tmp) <- c("V1", "NM", "chrom", "strand", "start", "end",
	##		   paste("V", 7:12, sep=""), "gene_name", paste("V", 14:16, sep=""))
	colClasses <- c("integer", "character", "character", "factor",
			"integer", "integer",
			"integer", "integer",
			"integer",
			"character", "character",
			"integer", rep("character", 4))
	tmp <- read.delim(filename, header=FALSE,
			  colClasses=colClasses)
	tmp <- tmp[, c(2:6, 13)]
	colnames(tmp) <- c("NM", "chrom", "strand", "start", "end", "gene_name")
	chrom <- sapply(tmp$chrom, function(x) strsplit(x, "chr")[[1]][2])
	tmp$chrom <- chromosome2integer(chrom)
	tmp <- tmp[!is.na(tmp$chrom), ]
	refGene <- RangedData(IRanges(tmp$start, tmp$end),
			      chrom=tmp$chrom,
			      strand=tmp$strand,
			      NM=tmp$NM,
			      gene_name=tmp$gene_name)
	refGene
}


filterCommonRegion <- function(deletion.ranges, CHR, FRAME=1e6){
	deletion.22 <- deletion.ranges[deletion.ranges$chrom == CHR, ]
	nn <- deletion.22$n.overlap
	index <- which(deletion.22$n.overlap == max(nn))
	##samples.22 <- strsplit(deletion.22$others[index], ", ")[[1]]
	samples.22 <- unique(as.character(sapply(deletion.22$others[index], function(x) strsplit(x, ", ")[[1]])))
	index.case <- deletion.22$id[index]
	samples.22 <- unique(c(index.case, samples.22))
	samples.22 <- samples.22[samples.22 != "NA"]
	deletion.22 <- deletion.22[deletion.22$id %in% samples.22, ]
##	deletion.22 <- deletion.22[!is.na(deletion.22$others), ]
	others <- unique(unlist(sapply(deletion.22$others, function(x) strsplit(x, ", ")[[1]])))
	others <- others[others != "NA"]
	##index <- which(deletion.ranges$
	deletion.22 <- deletion.22[deletion.22$id %in% others, ]
	##deletion.22 <- deletion.22[deletion.22$others != "NA", ]
	samples.22 <- unique(deletion.22$id)
	## if samples has just one overlap, see if its within 200kb of the sample with most overlap
	position.22 <- c(min(start(deletion.22)), max(end(deletion.22)))
	return(list(deletion.22, position.22))
}

plotCytobandWithRanges <- function(deletion.ranges, CHR, xlim, FRAME=1e6, ...){
	deletion.22 <- deletion.ranges[deletion.ranges$chrom == CHR, ]
	nn <- deletion.22$n.overlap
	index <- which(deletion.22$n.overlap == max(nn))
	##samples.22 <- strsplit(deletion.22$others[index], ", ")[[1]]
	samples.22 <- unique(as.character(sapply(deletion.22$others[index], function(x) strsplit(x, ", ")[[1]])))
	index.case <- deletion.22$id[index]
	samples.22 <- unique(c(index.case, samples.22))
	samples.22 <- samples.22[samples.22 != "NA"]
	deletion.22 <- deletion.22[deletion.22$id %in% samples.22, ]
##	deletion.22 <- deletion.22[!is.na(deletion.22$others), ]
	others <- unique(unlist(sapply(deletion.22$others, function(x) strsplit(x, ", ")[[1]])))
	others <- others[others != "NA"]
	##index <- which(deletion.ranges$
	deletion.22 <- deletion.22[deletion.22$id %in% others, ]
	##deletion.22 <- deletion.22[deletion.22$others != "NA", ]
	samples.22 <- unique(deletion.22$id)
	## if samples has just one overlap, see if its within 200kb of the sample with most overlap
	position.22 <- c(min(start(deletion.22)), max(end(deletion.22)))
	if(missing(xlim)){
		xlim <- c(min(position.22) -FRAME,
			  max(position.22) +FRAME)
	}
##	cyto.coords <-
##	plotCytoband(CHR, label.cytoband=FALSE, cytoband.ycoords=c(0, 0.1), xlim=xlim,
##		     ylim=c(0,length(samples.22)+0.5))
	plot(0:1, 0:1, xlim=xlim,ylim=c(0.5, length(samples.22)+0.5),
	     type="n", xlab="", ylab="", yaxt="n", xaxt="n")
	ii <- seq(0.15, 1, by=(1-0.15)/length(samples.22))
	h <- 0.3
	for(i in seq_along(samples.22)){
		this.range <- deletion.22[deletion.22$id == samples.22[i], ]
		## for each row of this subject, draw a polygon.
		for(j in 1:nrow(this.range)){
			x <- c(start(this.range)[j],end(this.range)[j])
			xx <- c(x, rev(x))
			y <- c(i-h, i-h, i+h, i+h)
			polygon(xx, y, col="grey60")
##			text(max(xx), mean(y), labels=paste("n =", this.range$num.mark[j]), col="blue")
			if(nrow(this.range) == 1)
				text(xlim[2], mean(y), labels=this.range$num.mark[j], col="grey30", cex=0.6, adj=1)
		}
		if(nrow(this.range) > 1){
			text(xlim[2], mean(y), labels=median(this.range$num.mark), adj=1, cex=0.6, col="grey30")
		}
	}
	axis(2, at=seq_along(samples.22), labels=samples.22, cex.axis=0.6, adj=0)
##	abline(v=position.22, col="blue", lwd=2, lty=2)
	text(0, mean(position.22), paste(diff(position.22)/1e3, "kb"))
	return(list(xlim=xlim, v=position.22))
}
