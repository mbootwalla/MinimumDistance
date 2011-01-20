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

doSegmentation <- function(){
	if(!exists("batch")) stop("batch variable should be defined in the submit_logR_cbs.R script")
	if(!exists("NN")) stop("batch size variable (NN) should be specified in submitter script")
	sample.index <- splitIndicesByLength(1:7599, NN)[[batch]]
	rD <- cbs(lrSet, sample.index=sample.index)
	save(rD, file=file.path(outdir, paste("rD_", batch, ".rda", sep="")))
	q("no")
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
				     space=rep(CHR, nrow(segmean_ranges)),
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


overlapsCentromere <- function(myranges, centromere.ranges, CHR){
	myranges.bak <- myranges
	centromere.ir <- IRanges(start(centromere.ranges)[CHR],
				 end(centromere.ranges)[CHR])
	if(!is(class(myranges), "IRanges")){
		myranges <- IRanges(start(myranges), end(myranges))
	}
	cnt <- countOverlaps(myranges, centromere.ir)
	return(cnt > 0)
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


minDistanceDeletion <- function(ranges, minDistanceSet, offspring.rule, CHR){
	mads <- minDistanceSet$Mad
	offspring.ids <- split(1:nrow(ranges), ranges$id)
	nn <- sapply(offspring.ids, length)
	thr <- rep(offspring.rule(mads), nn)
	is.deletion <- ifelse(ranges$seg.mean > thr, TRUE, FALSE)

	B <- baf(minDistanceSet)
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
	marker.index <- unique(unlist(marker.index))
	marker.index <- marker.index[order(marker.index)]
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
	stopifnot(identical(featureNames(minDistanceSet)[I], featureNames(bsSet)[I]))
	if(missing(id)){
		sample.names <- unique(substr(ranges1$id, 1, 5))
	} else sample.names <- substr(id, 1, 5)
	father.names <- paste(sample.names, "03", sep="_")
	mother.names <- paste(sample.names, "02", sep="_")
	father.index <- match(father.names, sampleNames(bsSet))
	mother.index <- match(mother.names, sampleNames(bsSet))
	offspr.index <- match(paste(sample.names, "01", sep="_"), sampleNames(bsSet))
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
	index <- match(ranges1$id, sampleNames(bsSet))
	mset$mad <- bsSet$MAD[index]
	rm(logR.F, logR.M, logR.O, baf.F, baf.M, baf.O, mindist); gc()
	return(mset)
}

collectAllRangesOfSize <- function(SIZE, bsSet, minDistanceSet,
				   minDistanceRanges,
				   outdir, MIN=1, MAX=4, lambda=0.1){
	library(SNPchip)
	data(chromosomeAnnotation)
	centromere.ranges <- GRanges(seqnames=Rle(paste("chr", 1:22, sep=""), rep(1,22)),
				     ranges=IRanges(chromosomeAnnotation[1:22, "centromereStart"],
				     chromosomeAnnotation[1:22, "centromereEnd"]))
	ranges2 <- minDistanceRanges[minDistanceRanges$seg.mean > 0 & minDistanceRanges$num.mark >= SIZE, ]
	index <- match(ranges2$id, sampleNames(minDistanceSet))
	open(minDistanceSet$MAD)
	mads <- minDistanceSet$MAD[index]
	x <- ranges2$num.mark
	p <- lambda*exp(-lambda*x)
	MIN <- 1; MAX <- 4
	b <- 1/(MAX - MIN)
	a <- MIN * b
	numberMads <- ((p-min(p))/(max(p)-min(p)) + a)/b
	thr <- numberMads * mads
	thr[thr < 0.2] <- 0.2
	ranges2$is.deletion <- ifelse(ranges2$seg.mean >= thr, TRUE, FALSE)
	deletion.RD <- ranges2[ranges2$is.deletion, ]
	deletion.RD$pHet <- NA
	deletion.RD$noCentromere.overlap <- NA
	deletion.ir <- IRanges(start(deletion.RD), end(deletion.RD))
	for(CHR in 1:22){
		cat("Chr ", CHR, "\n")
		i1 <- which(chromosome(bsSet) == CHR)
		i2 <- which(deletion.RD$chrom == CHR)

		fD <- fData(bsSet)[i1, ]
		fd.ir <- IRanges(fD$position-12, fD$position+12)

		##for each range, get the indices of probes for which it overlaps
		tmp <- matchMatrix(findOverlaps(deletion.ir[i2, ], fd.ir))

		## split the indices for the matching probes by the range
		index.list <- split(tmp[, "subject"], tmp[, "query"])
		fns.list <- lapply(index.list, function(i, fns) fns[i], fns=rownames(fD))

		marker.index <- i1[unique(as.integer(unlist(index.list)))]
		sample.index <- match(deletion.RD$id[i2], sampleNames(bsSet))
		B <- as.matrix(baf(bsSet)[marker.index, sample.index])
		pHet <- rep(NA, length(sample.index))
		for(i in seq_along(index.list)){
			ii <- match(fns.list[[i]], rownames(B))
			b <- B[ii, i]
			pHet[i] <- mean(b > 0.2 & b < 0.8, na.rm=TRUE)
		}
		deletion.RD$pHet[i2] <- pHet
		deletion.RD$noCentromere.overlap[i2] <- !(overlapsCentromere(deletion.RD[i2, ], centromere.ranges, CHR))
	}
	deletion.RD <- deletion.RD[order(deletion.RD$chrom, start(deletion.RD)), ]
	return(deletion.RD)
}

## for denovo-amplifications
collectAllRangesOfSize2 <- function(SIZE, bsSet,
				    maxDistanceSet,
				    maxDistanceRanges,
				    outdir, MIN=1, MAX=4,
				    lambda=0.1, upper.limit=-0.5){
	data(chromosomeAnnotation)
	centromere.ranges <- GRanges(seqnames=Rle(paste("chr", 1:22, sep=""), rep(1,22)),
				     ranges=IRanges(chromosomeAnnotation[1:22, "centromereStart"],
				     chromosomeAnnotation[1:22, "centromereEnd"]))
	##ranges2 <- maxDistanceRanges[maxDistanceRanges$seg.mean < 0 & maxDistanceRanges$num.mark >= SIZE, ]
	ranges2 <- maxDistanceRanges[maxDistanceRanges$seg.mean < upper.limit & maxDistanceRanges$num.mark >= SIZE, ]
	index <- match(ranges2$id, sampleNames(maxDistanceSet))
	open(maxDistanceSet$MAD)
	mads <- maxDistanceSet$MAD[index]
	x <- ranges2$num.mark
	p <- lambda*exp(-lambda*x)
	MIN <- 1; MAX <- 4
	b <- 1/(MAX - MIN)
	a <- MIN * b
	numberMads <- ((p-min(p))/(max(p)-min(p)) + a)/b
	thr <- -numberMads * mads
	thr[thr > upper.limit] <- upper.limit
	ranges2$is.altered <- ifelse(ranges2$seg.mean <= thr, TRUE, FALSE)
	altered.RD <- ranges2[ranges2$is.altered, ]
	altered.RD <- altered.RD[order(altered.RD$chrom, start(altered.RD)), ]
	return(altered.RD)
}

getRefGene <- function(filename="~/Data/Downloads/hg18_refGene.txt"){
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

getLocalId <- function(rd.object){
	idmap <- read.csv("~/projects2/Beaty/inst/extdata/gwas_to_local_id.csv", stringsAsFactors=FALSE)
	idmap <- idmap[idmap$cidr_name %in% rd.object$id, ]
	if(any(!rd.object$id %in% idmap$cidr_name)){
		print("The following ids in rd.object are not in the gwas_to_local_id file:")
		print(unique(rd.object$id[!rd.object$id %in% idmap$cidr_name]))
		rd.object <- rd.object[rd.object$id %in% idmap$cidr_name, ]
	}
	index <- match(rd.object$id, idmap$cidr_name)
	idmap <- idmap[index, ]
	stopifnot(identical(idmap$cidr_name, rd.object$id))
	rd.object$local_id <- idmap$local_id2
	return(rd.object)
}



findSubjectsInRange <- function(object, ##RangedData
				FRAME.SHIFT=200e3,
				MIN.SIZE=10){
	object <- object[order(object$num.mark, decreasing=T), ]
	object <- object[object$num.mark >= MIN.SIZE, ]
	object <- object[order(object$chrom), ]
	object$others <- NA
	chrom <- unique(object$chrom)
	for(i in seq_along(chrom)){
		CHR <- chrom[i]
		cat("Chr ", CHR, "\n")
		chr.range <- object[object$chrom == CHR, ]
		##segmean_ranges <- getSegMeans(outdir, CHR=CHR)
		load(file.path(beadstudiodir(), paste("cbs_chr" , CHR, ".rda", sep="")))
		ranges2 <- object[object$chrom == CHR, ]
		##Look for all ranges that have a denovo event in the region
		query.ir <- IRanges(start(chr.range) - FRAME.SHIFT,
				    end(chr.range) + FRAME.SHIFT)
		matching.subjects <- list()
		for(j in 1:nrow(chr.range)){
			## subjects -- all subjects except for the index case
			ranges3 <- ranges2[ranges2$id != chr.range$id[j], ]
			subj.ir <- IRanges(start(ranges3), end(ranges3))
			tmp <- matchMatrix(findOverlaps(query.ir[j, ], subj.ir))
			if(nrow(tmp) == 0) {
				matching.subjects[[j]] <- NA
			} else{
				matching.subjects[[j]] <- unique(ranges3$id[tmp[, 2]])
			}
		}
		tmp <- unlist(lapply(matching.subjects, function(x) paste(x, collapse=", ")))
		index <- which(object$chrom==CHR)
		object$others[index] <- tmp
	}
	matching.subjects <- object$others
	ms <- lapply(matching.subjects, function(x) strsplit(x, ", ")[[1]])
	##n.overlap <- sapply(matching.subjects, length)
	n.overlap <- sapply(ms, length)
	n.overlap[is.na(ms) | ms == "NA"] <- 0
	##object$others <- tmp
	object$n.overlap <- n.overlap
	return(object)
}

getDnaSource <- function(bsSet){
	dna <- as.character(bsSet$DNA.Source)
	cell.lines <- c(grep("CEPH", dna),
			grep("HAN CHINESE", dna),
			grep("JAPANESE", dna),
			grep("YORUBA", dna))
	dna[cell.lines] <- "HapMap"
	which.wga <- grep("WGA", dna)
	dna[which.wga] <- "WGA"
	dna[dna=="dried blood spot"] <- "blood"
	which.mouth <- match("mouth", dna)
	dna[which.mouth] <- "saliva"
	return(dna)
}
