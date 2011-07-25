##catFun <- function(rd.query, rd.subject, same.id=TRUE){
##	##browser()
##	##stopifnot(nrow(rd.query) == 1)
##	##rd.s <- rd.subject[seq(length=size), ]
####	if(!any(chromosome(rd.query) %in% chromosome(rd.subject))) {
####		return(0)
####	} else{
##	CHR <- chromosome(rd.query)
##	index <- which(chromosome(rd.query) == chromosome(rd.subject) & sampleNames(rd.query)==sampleNames(rd.subject)
##	rd.s <- rd.subject[chromosome(rd.subject) == CHR, ]
##
##
####	if(same.id){
##	index <- which(rd.subject$id == rd.query$id)
##	if(length(index) > 0){
##		rd.s <- rd.subject[index, ]
##	} ##else return(0)
####}
##	ir.q <- IRanges(start(rd.query),end(rd.query))
##	ir.s <- IRanges(start(rd.s),end(rd.s))
##	count <- countOverlaps(ir.q, ir.s)
##	return(count)
##}
catFun2 <- function(rd.query, rd.subject){
	##stopifnot(nrow(rd.query) == nrow(rd.subject)) ## must compare same list size
	ir.q <- IRanges(start(rd.query), end(rd.query))
	ir.s <- IRanges(start(rd.subject), end(rd.subject))
	mm <- findOverlaps(ir.q, ir.s)
	query.index <- queryHits(mm)
	subject.index <- subjectHits(mm)
	index <- which(chromosome(rd.query)[query.index] == chromosome(rd.subject)[subject.index] &
			sampleNames(rd.query)[query.index] == sampleNames(rd.subject)[subject.index])
	if(length(index) > 0){
		query.index <- unique(query.index[index])
		p <- length(query.index)/nrow(rd.query)
		if(p > 1) browser()
	} else {
		p <- 0
	}
	return(p)
}

discAtTop <- function(rd.query, rd.subject, verbose=TRUE){
	ir.q <- IRanges(start(rd.query), end(rd.query))
	ir.s <- IRanges(start(rd.subject), end(rd.subject))
	mm <- findOverlaps(ir.q, ir.s)
	query.index <- queryHits(mm)
	subject.index <- subjectHits(mm)
	index <- which(chromosome(rd.query)[query.index] == chromosome(rd.subject)[subject.index] &
		       sampleNames(rd.query)[query.index] == sampleNames(rd.subject)[subject.index])
	query.index <- unique(query.index[index])
	##subject.index <- unique(subject.index[index])
	notOverlapping.index <- seq(length=nrow(rd.query))[!seq(length=nrow(rd.query)) %in% query.index]
	rd.query[notOverlapping.index, ]
}

concAtTop <- function(ranges.query, ranges.subject, list.size, verbose=TRUE){
	p <- rep(NA, length(list.size))
	pAny1 <- rep(NA, length(list.size))
	pAny2 <- rep(NA, length(list.size))
	if(verbose) {
		message("Calculating the proportion of ranges in common for list sizes 1 to ", max(list.size))
		pb <- txtProgressBar(min=0, max=length(p), style=3)
	}
	for(i in seq_along(list.size)){
		if(verbose) setTxtProgressBar(pb, i)
		L <- list.size[i]
		p[i] <- catFun2(ranges.query[seq(length=L), ], ranges.subject[seq(length=L), ])
		pAny1[i] <- catFun2(ranges.query[seq(length=L), ], ranges.subject)
		pAny2[i] <- catFun2(ranges.subject[seq(length=L), ], ranges.query)
	}
	if(verbose) close(pb)
	return(list(p=p, pAny.queryList=pAny1, pAny.subjectList=pAny2))
}


##concAtTop <- function(ranges.query, ranges.subject, listSize, state,
##		      same.id=TRUE,
##		      verbose=TRUE, msg=".", method){
##	stopifnot(!missing(method))
##	stopifnot(length(method)==2)
##	if(!missing(state)){
##		if(verbose) message("Subsetting query and subject ranges by state ", state)
##		ranges.subject <- ranges.subject[state(ranges.subject) == state, ]
##		ranges.query <- ranges.query[state(ranges.query) == state, ]
##	} else{
##		if(verbose) message("State not specified. Checking concordance for denovo call")
##	}
##	stopifnot(listSize <= nrow(ranges.subject) & listSize <= nrow(ranges.query))
####	if(!same.id){
####		## identify the biggest hits in a region.  Carry only
####		## the biggest hit forward.
####		## -- can do this by matching to itself
####		message("Matching by locus...")
####		message("\tIdentifying unique query regions ")
####		nn <- coverage(ranges.query)
####		ranges.by.size <- ranges.query[order(nn, decreasing=TRUE), ]
####		##r1 <- ranges.by.size
####		NR <- nrow(ranges.by.size)
####		i <- 1
####		while(i < NR){
####			overlaps <- findOverlaps(ranges.by.size[i, ], ranges.by.size)
####			subj.index <- subjectHits(overlaps)[-1]
####			##make sure the chromosome is the same
####			subj.index <- subj.index[chromosome(ranges.by.size[i, ]) == chromosome(ranges.by.size)[subj.index]]
####			if(length(subj.index) > 0)
####				ranges.by.size <- ranges.by.size[-subj.index, ]
####			NR <- nrow(ranges.by.size)
####			i <- i+1
####		}
####		ranges.query <- ranges.by.size
####		message("\tIdentifying unique subject regions ")
####		nn <- coverage(ranges.subject)
####		ranges.by.size <- ranges.subject[order(nn, decreasing=TRUE), ]
####		NR <- nrow(ranges.by.size)
####		i <- 1
####		while(i < NR){
####			overlaps <- findOverlaps(ranges.by.size[i, ], ranges.by.size)
####			subj.index <- subjectHits(overlaps)[-1]
####			if(length(subj.index) > 0)
####				ranges.by.size <- ranges.by.size[-subj.index, ]
####			NR <- nrow(ranges.by.size)
####			i <- i+1
####		}
####		ranges.subject <- ranges.by.size
####	}
##	if(verbose) message("Ranking query and subject ranges by coverage")
##	ranges.subject$rank <- rank(-coverage(ranges.subject), ties.method="min")
##	ranges.query$rank <- rank(-coverage(ranges.query), ties.method="min")
##	listSize <- min(listSize, min(nrow(ranges.query), nrow(ranges.subject)))
##	top.query <- ranges.query[order(ranges.query$rank, decreasing=FALSE)[1:listSize], ]
##
##	countAny <- rep(NA, listSize)
##	##rankInSubject <- rep(NA, nrow(top.query))
##	if(verbose) {
##		message("Calculating the proportion of ranges in common for list sizes 1 to ", listSize)
##		pb <- txtProgressBar(min=0, max=length(countAny), style=3)
##	}
##	for(i in seq(length=nrow(top.query))){
##		if(verbose) setTxtProgressBar(pb, i)
##		if(i %% 100 == 0) cat(".")
##		countAny[i] <- catFun(rd.query=top.query[i, ], rd.subject=ranges.subject, same.id=same.id)
##	}
##	Iany <- countAny > 0
##	message("returning proportion in common (p) and coverage in query (cov)")
##	pAny <- sapply(1:nrow(top.query), function(x, I) mean(I[1:x]), I=Iany)
##
##	top.subject <- ranges.subject[order(ranges.subject$rank, decreasing=FALSE)[1:listSize], ]
##	count <- rep(NA, nrow(top.query))
##	count2 <- rep(NA, nrow(top.query))
##	##rankInSubject <- rep(NA, nrow(top.query))
##	if(verbose) {
##		message("Calculating the proportion of ranges in common for list sizes 1 to ", listSize)
##		pb <- txtProgressBar(min=0, max=length(count), style=3)
##	}
##	for(i in seq(length=nrow(top.query))){
##		if(verbose) setTxtProgressBar(pb, i)
##		if(i %% 100 == 0) cat(".")
##		count[i] <- catFun(rd.query=top.query[seq(length=i), ], rd.subject=top.subject[seq(length=i), ], same.id=same.id)
##		count2[i] <- catFun(rd.query=top.subject[seq(length=i), ], rd.subject=top.query[seq(length=i), ], same.id=same.id)
##		if(count[i] != count2[i]) browser()
##	}
##	## which ranges in query were not found in subject
##	ir1 <- IRanges(start(top.query), end(top.query))
##	chr1 <- chromosome(top.query)
##	id1 <- sampleNames(top.query)
##	ir2 <- IRanges(start(top.subject), end(top.subject))
##	chr2 <- chromosome(top.subject)
##	id2 <- sampleNames(top.subject)
##	mm <- findOverlaps(ir1, ir2)
##	##index <- which(chr1==chr2 & id1 == id2)
##	query.hits <- queryHits(mm)
##	subject.hits <- subjectHits(mm)
##	index <- which(chr1[query.hits] == chr2[subject.hits] & id1[query.hits] == id2[subject.hits])
##	query.hits <- query.hits[index]
##	subject.hits=subject.hits[index]
##	## ranges in both
##	top.query$method <- method[1]
##	top.subject$method <- method[2]
##	if(length(query.hits) > 1){
##		concordant=list(query=top.query[query.hits,],
##		subject=top.subject[subject.hits,])
##	}
##	index2 <- seq(length=nrow(top.query))[-query.hits]
##	if(length(index2) > 1){
##		discordant=list(queryNotInSubject=top.query[index2, ],
##		subjectNotInQuery=top.subject[index2, ])
##	}
##	if(verbose) close(pb)
##	I <- count > 0
##	message("returning proportion in common (p) and coverage in query (cov)")
##	p <- sapply(1:nrow(top.query), function(x, I) mean(I[1:x]), I=I)
##	cov <- coverage(top.query)
##	return(list(pAny=pAny,
##		    p=p, cov=cov,
##		    concordant=concordant,
##		    discordant=discordant))
##}

##assessConcordance <- function(md332, penn332, listSize=500){
##	ct332 <-  concAtTop(md332, penn332, listSize=listSize, state="332", same.id=FALSE)
##	## minimum distance | pennCNV
##	ct332r <- concAtTop(penn332, md332, listSize=listSize, state="332", same.id=FALSE)
##	ct332.2 <-  concAtTop(md332, penn332, listSize=listSize, state="332")
##	ct332r.2 <- concAtTop(penn332, md332, listSize=listSize, state="332")
##	ct <- list(ct332, ct332r, ct332.2, ct332r.2)
##	return(ct)
##}

notCalled <- function(ranges.query, ranges.subject, listSize, sample.match=TRUE){
	message("Returning ranges in query that do not overlap with top ranges in subject")
	ranges.subject$rank <- rank(-coverage(ranges.subject), ties.method="min")
	ranges.query$rank <- rank(-coverage(ranges.query), ties.method="min")
	message("Assessing top ", listSize, " query ranges for hit in subject")
	top.query <- ranges.query[order(ranges.query$rank, decreasing=FALSE)[1:listSize], ]
	ranges.subject <- ranges.subject[order(ranges.subject$rank,decreasing=FALSE)[1:listSize], ]
	##message("Looking at all subject ranges regardless of coverage")
	##top.subject <- ranges.subject[order(ranges.subject$rank, decreasing=FALSE)[1:listSize], ]
	##top.subject <- ranges.subject[order(ranges.subject$rank, decreasing=FALSE)[1:listSize], ]
	overlap <- findOverlaps(top.query, ranges.subject)
	subj.index <- subjectHits(overlap)
	quer.index <- queryHits(overlap)
	## what are the chromosomes for the subject hits
	chrom.subj <- ranges.subject$chrom[subj.index]
	chrom.quer <- top.query$chrom[quer.index]
	id.subj <- ranges.subject$id[subj.index]
	id.quer <- top.query$id[quer.index]
	if(sample.match){
		##	## eliminate those for which the chromosome is not the same
		message("Requiring sample id to match")
		ii <- which(chrom.subj == chrom.quer & id.subj==id.quer)
	} else {
		message("Requiring only overlap of the region on the same chromosome")
		ii <- which(chrom.subj == chrom.quer)
	}
	subj.index <- subj.index[ii]
	## index of ranges in query that have a match
	quer.index <- quer.index[ii]
	absent.index <- seq(length=nrow(top.query))[!seq(length=nrow(top.query)) %in% quer.index]
	return(top.query[absent.index, ])
}

## from Smyth,2004
estimateDnot <- function(s2.g, d.g, epsilon=1e-8){
	## d.g can be allowed to differ between genes
	## doesn't d.g = n-1?
	stopifnot(length(s2.g) > 1) ## vector
	z.g <- log(s2.g)
	e.g <- z.g - digamma(d.g/2) + log(d.g/2)
	e.bar <- mean(e.g, na.rm=TRUE)
	x <- mean((e.g-e.bar)^2*(d.g+1)/d.g - trigamma(d.g/2), na.rm=TRUE)
	## solve newton raphson
	y.previous <- 0.5 + 1/x
	tetragamma <- function(x) psigamma(x, 2)
	## digamma(x)=psigamma(x, 0)
	## trigamma(x)=psigamma(x, 1)
	## tetragamma(x)=psigamma(x, 2)
	err <- 1
	counter <- 1
	while(err > epsilon){
		delta <- trigamma(y.previous)*(1-trigamma(y.previous)/x)/tetragamma(y.previous)
		y.next <- y.previous+delta
		err <- -delta/y.next
		##cat(err, "\n")
		counter <- counter+1
		if(counter > 100) {
			warning("Newton-Raphson for d.0 did not converge in 100 iterations")
			break()
		}
		y.previous <- y.next
	}
	d.0 <- y.previous*2
	s2.0 <- exp(e.bar + digamma(d.0/2) - log(d.0/2))
	res <- c(d.0, s2.0)
	names(res) <- c("d.0", "s2.0")
	return(res)
}

correspondingCall <- function(ranges.query, ranges.subject, subject.method){
	overlap <- findOverlaps(ranges.query, ranges.subject)
	subj.index <- subjectHits(overlap)
	quer.index <- queryHits(overlap)
	## what are the chromosomes for the subject hits
	index <- which(chromosome(ranges.query)[quer.index] == chromosome(ranges.subject)[subj.index] &
		       sampleNames(ranges.query)[quer.index] == sampleNames(ranges.subject)[subj.index])
	if(length(index) == 0) return("no overlap")
	matching.index <- subj.index[index]
	res <- ranges.subject[matching.index, ]
##	chrom.subj <- ranges.subject$chrom[subj.index]
##	chrom.quer <- ranges.query$chrom[quer.index]
##	id.subj <- ranges.subject$id[subj.index]
##	id.quer <- ranges.query$id[quer.index]
##	##	## eliminate those for which the chromosome is not the same
##	ii <- which(chrom.subj == chrom.quer & id.subj==id.quer)
##	subj.index <- subj.index[ii]
##	## index of ranges in query that have a match
##	res <- ranges.subject[subj.index, ]
	if(!missing(subject.method)) res$method <- subject.method
	return(res)
}






constructFeatureData <- function(file, cdfName){
	dat <- read.delim(file, colClasses=c("character", "numeric", "numeric"))
	mito.index <- grep("Mito", dat$Name)
	dat <- dat[-mito.index, ]
 	featureData <- oligoClasses:::featureDataFrom(paste(cdfName, "Crlmm", sep=""))
	stopifnot(nrow(featureData) == nrow(dat))
	return(featureData)
}

constructLogRatioSet <- function(fnames, filename, cdfName, samplesheet){
	featureData <- constructFeatureData(fnames[1], cdfName)
	nc <- length(fnames)
	nr <- nrow(dat)
	baf <- initializeBigMatrix("baf", nr, nc, vmode="double")
	logR <- initializeBigMatrix("logR", nr, nc, vmode="double")
	setMethod("annotatedDataFrameFrom", "ffdf", Biobase:::annotatedDataFrameFromMatrix)
	bsSet <- new("LogRatioSet", logRRatio=logR, BAF=baf, annotation=cdfName)
	position.order <- order(featureData$chromosome, featureData$position)
	featureData <- featureData[position.order, ]
	featureData(bsSet) <- featureData
	featureNames(bsSet) <- featureNames(featureData)
	sns <- basename(fnames)##, 1,8)
	sns.short <- substr(sns, 1, 8)
	samplesheet$short.name <- substr(samplesheet$Sample.Name, 1, 8)
	index <- match(sns.short, samplesheet$short.name)
	samplesheet <- samplesheet[index, ]
	stopifnot(all.equal(samplesheet$short.name, sns.short))

	dup.index <- which(duplicated(sns))
	sns[dup.index] <- paste(sns[dup.index], 1:length(dup.index), sep="_")
	sampleNames(bsSet) <- sns
	##samplesheet <- samplesheet[match(sampleNames(bsSet), samplesheet$Sample.Name), ]
	##stopifnot(identical(samplesheet$Sample.Name, sampleNames(bsSet)))
	pData(bsSet) <- samplesheet
	sampleNames(phenoData(bsSet)) <- sns
	mads <- initializeBigVector(name="mads", n=ncol(bsSet), vmode="double")
	bsSet$MAD <- mads
}

getChromosomeArm <- function(chrom, pos){
	if(!is.integer(chrom)) {
		chrom <- chromosome2integer(chrom)
	}
	if(!all(chrom %in% 1:24)){
			warning("Chromosome annotation is currently available for chromosomes 1-22, X and Y")
			##message("Please add/modify data(chromosomeAnnotation, package='SNPchip') to accomodate special chromosomes")
			pos <- pos[chrom%in%1:24]
			chrom <- chrom[chrom%in%1:24]
	}
	data(chromosomeAnnotation, package="SNPchip", envir=environment())
	chromosomeAnnotation <- as.matrix(chromosomeAnnotation)
	chrAnn <- chromosomeAnnotation
	uchrom <- unique(SNPchip:::integer2chromosome(chrom))
	chromosomeArm <- vector("list", length(uchrom))
	positionList <- split(pos, chrom)
	positionList <- positionList[match(unique(chrom), names(positionList))]
	for(i in seq_along(unique(chrom))){
		chromosomeArm[[i]] <- ifelse(positionList[[i]] <= chrAnn[uchrom[i], "centromereEnd"], "p", "q")
	}
	chromosomeArm <- unlist(chromosomeArm)
	chrom <- paste(chrom, chromosomeArm,sep="")
}


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

loadRangesCbs <- function(outdir, pattern, CHR, name){
	fname <- list.files(outdir, pattern=pattern, full.name=TRUE)
	if(missing(name)) stop("must specify object name")
	stopifnot(length(fname) == 1)
	load(fname)
	cbs.segs <- get(name)
	cbs.segs
}

## load and try to rbind the list
loadBatchFiles <- function(outdir, pattern, name){
	fnames <- list.files(outdir, pattern=pattern, full.name=TRUE)
	dfl <- vector("list", length(fnames))
	for(i in seq_along(fnames)){
		load(fnames[i])
		dfl[[i]] <- get(name)
	}
	if(length(dfl) == 1) dfl <- dfl[[1]]
	##if(length(dfl) > 1) df <- tryCatch(df <- do.call("rbind", dfl), error=function(e) NULL)
	return(dfl)
}

featuresInRange <- function(object, range, ...){
	featuresInXlim(object, start=start(range), end=end(range), CHR=range$chrom, ...)
}

addIndicesFromFeatureData <- function(rd.object, fD, FRAME=0){
	uchrom.rd <- unique(rd.object$chrom)
	uchrom.fd <- unique(fD$chromosome)
	stopifnot(length(uchrom.rd)==1)
	stopifnot(length(uchrom.fd)==1)
	stopifnot(uchrom.rd == uchrom.fd)
	pos.chr <- fD$position
	subj <- IRanges(pos.chr-25, pos.chr+25)
	tree <- IntervalTree(subj)
	query <- IRanges(start(rd.object)-FRAME, end(rd.object)+FRAME)
	tmp <- findOverlaps(query, tree)
	mm <- matchMatrix(tmp)
	indexlist <- split(mm[, 2], mm[, 1])
	stopifnot(length(indexlist) == nrow(rd.object))
	tmp <- lapply(indexlist, range)
	tmp <- do.call(rbind, tmp)
	rd.object$start.index <- tmp[,1]
	rd.object$end.index <- tmp[,2]
	rd.object
}

addIndicesForRanges <- function(rd.object, bsSet, FRAME=0, xlim,
				index.in.chromosome=FALSE, verbose=TRUE){
	if(verbose) {
		if(index.in.chromosome){
			message("adding start.index and end.index with respect to chromosome")
		} else {
			message("adding start.index and end.index with respect to index in bsSet object")
		}
	}
	stopifnot(all(order(chromosome(bsSet), position(bsSet))) > 0)
	pos <- position(bsSet)
	chrom <- chromosome(bsSet)
	rd.object$start.index <- NA
	rd.object$end.index <- NA
	J <- unique(rd.object$chrom)
	for(CHR in J){
		##space.index <- which(space(rd.object)==CHR)
		space.index <- which(rd.object$chrom==CHR)
		ranges.chr <- rd.object[space.index, ]
		chr.index <- which(chrom == CHR)
		pos.chr <- pos[chr.index]
		subj <- IRanges(pos.chr-12, pos.chr+12)
		tree <- IntervalTree(subj)
		if(missing(xlim)){
			query <- IRanges(start(ranges.chr)-FRAME, end(ranges.chr)+FRAME)
		} else query <- IRanges(xlim[1], xlim[2])
		tmp <- findOverlaps(query, tree)
		mm <- matchMatrix(tmp)
		## these are the indices for a specific chromosome.  Need to intersect
		## with chr.index to map back to the ori
		if(!index.in.chromosome){
			indicesList <- split(chr.index[mm[, 2]], mm[, 1])
			indexList <- lapply(indicesList, range)
		} else {
			indicesList <- split(mm[, 2], mm[, 1])
			indexList <- lapply(indicesList, range)
		}
		index <- do.call("rbind", indexList)
		rd.object$start.index[space.index] <- index[,1]
		rd.object$end.index[space.index] <- index[, 2]
	}
	rd.object
}

featuresInXlim <- function(object, start, end, CHR, pos, chrom, FRAME=0, FRAME.LEFT, FRAME.RIGHT){
	stopifnot(!missing(CHR))
	if(missing(FRAME.LEFT)) FRAME.LEFT <- FRAME
	if(missing(FRAME.RIGHT)) FRAME.RIGHT <- FRAME
	if(missing(start)) start <- 0
	if(missing(end)){
		require(SNPchip)
		data(chromosomeAnnotation)
		end <- chromosomeAnnotation[CHR, "chromosomeSize"]
	}
	start <- start-FRAME.LEFT
	end <- end+FRAME.RIGHT
	if(missing(pos) | missing(chrom)){
		which(position(object) >= start & position(object) <= end & chromosome(object) == CHR)
	} else which(pos >= start & pos <= end & chrom == CHR)
}

framePositionIndex <- function(object,  ##LogRatioSet or something similar
			       FRAME){  ##basepairs
		min.pos <- min(pos)-WINDOW.SIZE
		max.pos <- max(pos)+WINDOW.SIZE
}

overlapsCentromere <- function(myranges){
	require(SNPchip)
	data(chromosomeAnnotation)
	centromere.ranges <- RangedData(IRanges(chromosomeAnnotation[, "centromereStart"],
						chromosomeAnnotation[, "centromereEnd"]),
					chrom=rownames(chromosomeAnnotation))
	myranges.bak <- myranges
	chrom <- unique(myranges$chrom)
	overlaps.centromere <- rep(NA, nrow(myranges))
	for(CHR in chrom){
		centromere.ir <- IRanges(start(centromere.ranges)[CHR],
					 end(centromere.ranges)[CHR])
		ix <- which(myranges$chrom==CHR)
		ir <- IRanges(start(myranges)[ix],
			      end(myranges)[ix])
		overlaps.centromere[ix] <- countOverlaps(ir, centromere.ir) > 0
	}
	return(overlaps.centromere)
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

combineRanges <- function(deletion.ranges, amp.ranges){
	state <- deletion.ranges$state
	hemizygous.states <- c("332", "432", "342")
	homozygous.states <- c("331", "321", "231", "431", "341", "441", "221")
	deletion.ranges <- deletion.ranges[state %in% hemizygous.states | state %in% homozygous.states, ]
	amp.ranges <- amp.ranges[, colnames(amp.ranges) %in% colnames(deletion.ranges)]
	index <- match(colnames(amp.ranges), colnames(deletion.ranges))
	deletion.ranges2 <- deletion.ranges[,  index]
	stopifnot(all.equal(colnames(deletion.ranges2), colnames(amp.ranges)))
	ranges.all <- RangedData(IRanges(c(start(deletion.ranges2), start(amp.ranges)),
					 c(end(deletion.ranges2), end(amp.ranges))),
				 id=c(deletion.ranges2$id, amp.ranges$id),
				 chrom=c(deletion.ranges2$chrom, amp.ranges$chrom),
				 num.mark=c(deletion.ranges2$num.mark, amp.ranges$num.mark),
				 seg.mean=c(deletion.ranges2$seg.mean, amp.ranges$seg.mean),
				 state=c(deletion.ranges2$state, amp.ranges$state))
	ranges.all
}


pruneByFactor <- function(range.object, f, verbose){
	rd <- list()
	id.chr <- paste(range.object$id, chromosome(range.object), sep="_")
	ff <- unique(id.chr)
	##for(i in seq_along(unique(range.object$id))){
	if(verbose){
		message("Pruning ", length(ff), " files.")
		pb <- txtProgressBar(min=0, max=length(ff), style=3)
	}
	for(i in seq_along(ff)){
		if(verbose) setTxtProgressBar(pb, i)
		##id <- unique(range.object$id)[i]
		##(index <- which(range.object$id == id))
		index <- which(id.chr==ff[i])
		##trace(combineRangesByFactor, browser)
		rd[[i]] <- combineRangesByFactor(range.object[index, ], f=f[index])
	}
	if(verbose) close(pb)
	ok <- tryCatch(tmp <- do.call("rbind", rd), error=function(e) FALSE)
	if(is(ok, "logical")) tmp <- rd
	return(tmp)
}

combineRangesByFactor <- function(range.object, f){
	i <- which(is.na(f))
	j <- 1
	while(length(i) > 0){
		if(is.na(f[1])){
			f[1] <- f[2]
		} else {
			f[is.na(f)] <- f[i-1]
		}
		i <- which(is.na(f))
		j <- j+1
		if(j > 10) stop("too many na's in f")
	}
	##stopifnot(all(!is.na(f)))
	ff <- cumsum(c(0, abs(diff(as.integer(as.factor(f))))))
	if(!any(duplicated(ff))) return(range.object)
	for(i in seq_along(unique(ff))){
		x <- unique(ff)[i]
		if(sum(ff==x) == 1) next()
		index <- which(ff==x)
		min.index <- min(index)
		max.index <- max(index)
		end(range.object)[index] <- max(end(range.object)[index])
		range.object$lik.state[index] <- sum(range.object$lik.state[index])
		range.object$end.index[index] <- max(range.object$end.index[index])
		range.object$seg.mean[index] <- sum((range.object$num.mark[index] * range.object$seg.mean[index]))/sum(range.object$num.mark[index])
		range.object$num.mark[index] <- sum(range.object$num.mark[index])
		j <- seq(length=nrow(range.object))
		index <- index[-1]
		j <- j[-index]
		if(length(j) == 0){
			stop()
		}
		ff <- ff[j]
		range.object <- range.object[j, ]
	}
	return(range.object)
}

madVsCoverage <- function(lambda=0.1, MIN=1, MAX=4, coverage=3:100){
	p <- lambda*exp(-lambda*coverage) ## 0 - 0.04 (Pr (X=x)
	b <- 1/(MAX - MIN)
	a <- MIN * b
	numberMads <- ((p-min(p))/(max(p)-min(p)) + a)/b
	list(x=coverage, y=numberMads)
}

getBedFiles <- function(xlim, CHR, rf, cnv, envir){
	objectname <- paste("rf.chr", CHR, sep="")
	needToLoad <- !crlmm:::isLoaded(objectname, envir)
	chr.name <- paste("chr", CHR, sep="")
	if(needToLoad){
		rf.chr <- rf[rf$chrom==chr.name, ]
		assign(objectname, rf.chr, envir=envir)
	}
	rf.chr <- crlmm:::getVarInEnv(objectname, envir)
	objectname <- paste("cnv.chr", CHR, sep="")
	needToLoad <- !crlmm:::isLoaded(objectname, envir)
	if(needToLoad){
		cnv.chr <- cnv[cnv$chrom==chr.name, ]
		assign(objectname, cnv.chr, envir=envir)
	}
	cnv.chr <- crlmm:::getVarInEnv(objectname, envir)
	rf.ir <- IRanges(rf.chr$txStart, rf.chr$txEnd)
	cnv.ir <- IRanges(cnv.chr$txStart, cnv.chr$txEnd)
	subj <- IRanges(xlim[1]*1e6, xlim[2]*1e6)
	rf.index <- matchMatrix(findOverlaps(rf.ir, subj))[, 1]
	cnv.index <- matchMatrix(findOverlaps(cnv.ir, subj))[, 1]
	rf.chr <- rf.chr[rf.index, ]
	cnv.chr <- cnv.chr[cnv.index, ]
	flatBed.genes <- flatten.bed(rf.chr)
	flatBed.genes$start <- flatBed.genes$start/1e3
	flatBed.genes$stop <- flatBed.genes$stop/1e3
	flatBed.cnv <- flatten.bed(cnv.chr)
	flatBed.cnv$start <- flatBed.cnv$start/1e3
	flatBed.cnv$stop <- flatBed.cnv$stop/1e3
	list(genes=flatBed.genes,
	     cnv=flatBed.cnv)
}

thresholdSegMeans <- function(ranges.object, ylim){
	ranges.object$seg.mean[ranges.object$seg.mean < ylim[1]] <- ylim[1]
	ranges.object$seg.mean[ranges.object$seg.mean > ylim[2]] <- ylim[2]
	ranges.object
}

combine.data.frames <- function(dist.df, penn.df){
	if(is.null(dist.df) & is.null(penn.df)) return(NULL)
	if(is.null(dist.df)) dist.df <- penn.df[integer(0), ]
	if(is.null(penn.df)) penn.df <- dist.df[integer(0), ]
	combined.df <- rbind(dist.df, penn.df)
	combined.df <- combined.df[order(combined.df$chr), ]
	return(combined.df)
}

deletionStates <- function(){
	st1 <- offspring.hemizygous()
	st2 <- offspring.homozygous()
	as.integer(c(st1,st2))
}
offspring.hemizygous <- function() c("332", "432", "342", "442")
offspring.homozygous <- function() c("331", "321", "231", "431", "341", "441", "221", "421")
duplicationStates <- function() as.integer(c("335", "334", "224", "225", "115", "114", "124", "125", "214", "215", "324", "325", "234", "235", "124", "125", "214", "215", "314", "315", "134", "135"))
duplicationStatesPenn <- function() as.integer(c("335", "225", "115", "125", "215", "325", "235", "125", "215", "315", "135"))
isDenovo <- function(states) states %in% c(duplicationStates(), deletionStates())

cbsSegsForRange <- function(ranges.object, ylim, envir){
	stopifnot(nrow(ranges.object) == 1)
	id <- substr(ranges.object$id[[1]], 1, 5)
	CHR <- ranges.object$chrom[[1]]
	cbs.segs <- cbsSegsForChrom(CHR, envir)
	segs.index <- which(substr(cbs.segs$id, 1, 5) %in% id)
	cbs.segs <- cbs.segs[segs.index, ]
	cbs.segs$seg.mean[cbs.segs$seg.mean < ylim[1]] <- ylim[1]
	cbs.segs$seg.mean[cbs.segs$seg.mean > ylim[2]] <- ylim[2]
	cbs.segs
}

cbsSegsForChrom <- function(CHR, envir){
	segname <- paste("cbs_segs", CHR , sep="")
	cbs.segs <- cbsLoader(segname, CHR, envir)
	cbs.segs
}


insertCentromereBreak <- function(range.object, insertWhere=c("centromereStart", "centromereEnd")){
	data(chromosomeAnnotation)
	stopifnot(insertWhere %in% c("centromereStart", "centromereEnd"))
	chrom <- unique(range.object$chrom)
	stopifnot(all(is.integer(chrom)))
	j <- which(colnames(chromosomeAnnotation) == insertWhere)
	if(insertWhere=="centromereStart"){
		index <- which(start(range.object) <= chromosomeAnnotation[chrom, 1] & end(range.object) >= chromosomeAnnotation[chrom, 1])
		stopifnot(length(index) == 1)
		##range.object <- range.object[index, ]
		end(range.object)[index] <- chromosomeAnnotation[chrom, 1]
	}
	if(insertWhere=="centromereEnd"){
		index <- which(start(range.object) <= chromosomeAnnotation[chrom, 2] & end(range.object) >= chromosomeAnnotation[chrom, 2])
		stopifnot(length(index) == 1)
		start(range.object)[index] <- chromosomeAnnotation[chrom, 2]
	}
	range.object[index, ]
}

addCentromereRanges <- function(centromere.ranges){
	range1 <- insertCentromereBreak(centromere.ranges, insertWhere="centromereStart")
	range2 <- insertCentromereBreak(centromere.ranges, insertWhere="centromereEnd")
	ranges <- c(range1, range2)
}

replaceCentromereRanges <- function(range.object){
	overlaps.centromere <- overlapsCentromere(range.object)
	index <- which(overlaps.centromere)
	if(length(index) == 0) {
		message("no ranges overlap centromere.  do nothing")
	} else{
		message(length(index), " ranges overlap centromere.  Replacing with ranges that do not span the centromere")
		ranges.updated <- addCentromereRanges(range.object[index, ])
		range.object <- range.object[-index, ]
	}
	range.object
}

musInRange <- function(query, cbs.segs, id, chr){
	index <- which(cbs.segs$id == id)
	cbs.sub <- cbs.segs[index, ]
	subj <- IRanges(start(cbs.sub), end(cbs.sub))
	mm <- matchMatrix(findOverlaps(query, subj))[,2]
	cbs.sub[mm, ]
}

##short
ssampleNames <- function(object) substr(sampleNames(object), 1, 8)
##short short
sssampleNames <- function(object) substr(sampleNames(object), 1, 5)
##short short names
ssnames <- function(x) ss(names(x))
## short names
snames <- function(x) s(names(x))
ss <- function(x) substr(x, 1, 5)
s <- function(x) substr(x, 1, 8)
calculateChangeSd <- function(coverage=1:500, lambda=0.05, a=0.2, b=0.025)
	a + lambda*exp(-lambda*coverage)/b

pruneMD <- function(genomdat,
		  range.object,
		  physical.pos,
		  ##trimmed.SD, ##
		  lambda=0.05,
		  MIN.CHANGE=0.1,
		  SCALE.EXP=0.02,
		  MIN.COVERAGE=3,
		  weighted=FALSE,
		  weights=NULL) {
	stopifnot(length(unique(range.object$id)) == 1)
	stopifnot(length(unique(range.object$chrom)) == 1)
	##change.SD <- trimmed.SD*change.SD
	genomdat <- as.numeric(genomdat)
	coverage <- range.object$num.mark
	trimmed.SD <- unique(range.object$mindist.mad)
	stopifnot(length(trimmed.SD)==1)
	coverage <- coverage[-length(coverage)]
	if(FALSE){
		numberSds <- calculateChangeSd(coverage=3:100, lambda=lambda, a=MIN.CHANGE, b=SCALE.EXP)
		graphics:::plot(3:100, y, ylab="number of MADs", xlab="coverage")
	}
	##thrSD <- calculateChangeSd(coverage, lambda, trimmed.SD, change.SD)
		##change.SD <- change.SD  ##Thresholds for right cutpoint
	##cpt.loc <- cumsum(lseg) ## indices of the cutpoints same as coverage.
	cpt.loc <- range.object$end.index
	sdundo <- TRUE
	while(sdundo) {
		k <- length(cpt.loc)
		if (k>1) {
			coverage <- diff(c(0, cpt.loc))
			coverage <- coverage[-length(coverage)]
			##
			##  number of sds as a function of coverage
			##  -- segments with high coverage have small y
			##
			requiredNumberSd <- calculateChangeSd(coverage=coverage, lambda=lambda, a=MIN.CHANGE, b=SCALE.EXP)
			##
			## number of standard deviations
			segments0 <- cbind(c(1,1+cpt.loc[-k]),cpt.loc)
			## median copy number for each segment
			segmed <- apply(segments0, 1, function(i,x) {median(x[i[1]:i[2]], na.rm=T)}, genomdat)
			## absolute copy number difference of adjacent segments
 			##adsegmed <- abs(diff(segmed))
			adsegmed <- abs(diff(segmed))
			## number of standard deviations of observed shift
			empiricalNumberSd <- adsegmed/trimmed.SD
			if(any(empiricalNumberSd < requiredNumberSd | coverage < MIN.COVERAGE)){
				## drop order: coverage then distance
				##i <- which(adsegmed < thrSD | coverage < MIN.COVERAGE)
				i <- which(empiricalNumberSd < requiredNumberSd | coverage < MIN.COVERAGE)
				if(length(i) > 1){
					i <- i[order(coverage[i], adsegmed[i], decreasing=FALSE)[1]]
				}
				cpt.loc <- cpt.loc[-i]
			} else {
				sdundo <- FALSE
			}
		} else {
			sdundo <- FALSE
		}
	}
	lseg <- diff(c(0,cpt.loc)) ## back to coverage
	## update segment means
	segmeans <- 0*lseg
	ll <- uu <- 0
	for(i in 1:length(lseg)) {
		uu <- uu + lseg[i]
		if (weighted) {
			segmeans[i] <- sum(genomdat[(ll+1):uu]*weights[(ll+1):uu])/sum(weights[(ll+1):uu])
		} else {
			segmeans[i] <- mean(genomdat[(ll+1):uu], na.rm=TRUE)
		}
		ll <- uu
	}
	segments0 <- cbind(c(1,1+cpt.loc[-k]),cpt.loc)
	starts <- physical.pos[segments0[, 1]]
	ends <- physical.pos[segments0[, 2]]
	id <- unique(range.object$id)
	RangedData(IRanges(starts, ends),
		   id=id,
		   chrom=unique(range.object$chrom),
		   num.mark=lseg,
		   seg.mean=segmeans,
		   start.index=segments0[,1],
		   end.index=segments0[,2],
		   mindist.mad=range.object$mindist.mad[1])
}
initializeBigArray <- function(name, dim, vmode="integer", initdata=NA){
	if(isPackageLoaded("ff")){
		if(prod(dim) > 2^31){
			##Need multiple matrices
			## -- use ffdf
			## How many samples per ff object
			S <- floor(2^31/nr - 1)
			## How many ff objects
			L <- ceiling(nc/S)
			name <- paste(name, 1:L, sep="_")
			resultsff <- vector("list", L)
			for(i in 1:(L-1)){  ## the Lth object may have fewer than nc columns
				resultsff[[i]] <- createFF(name=name[i],
							   dim=c(nr, S),
							   vmode=vmode, initdata=initdata)
			}
			##the Lth element
			leftOver <- nc - ((L-1)*S)
			resultsff[[L]] <- createFF(name=name[L],
						   dim=c(nr, leftOver),
						   vmode=vmode, initdata=initdata)
			results <- do.call(ffdf, resultsff)
			rm(resultsff); gc()
		} else {
			results <- createFF(name=name,
					    dim=c(nr, nc),
					    vmode=vmode, initdata=initdata)
		}
	}  else {
		init <- switch(vmode,
			       integer=as.integer(initdata),
			       double=as.double(initdata),
			       character=as.character(initdata),
			       stop("Mode ", vmode, " not implemented for regular matrices"))
		results <- matrix(init, nr, nc)
	}
	return(results)
}

## pdf of standard normal
## the msm package has this stuff, but it seemed slow...
phi <- function(x, mu, sigma) dnorm(x, mu, sigma)
## cdf of standard normal
Phi <- function(x, mu, sigma) pnorm(x, mu, sigma)
## pdf of truncated normal on support [0, 1]
##tnorm <- function(x, mu, sigma) phi(x, mu, sigma)/(Phi(1, mu, sigma)-Phi(0, mu, sigma))
tnorm <- function(x, mean, sd, lower=0, upper=1){
	res <- phi(x, mean, sd)/(Phi(upper, mean, sd)-Phi(lower, mean, sd))
	ind <- which(x < lower | x > upper)
	if(any(ind)){
		res[ind] <- 0
	}
	res
}
TN <- tnorm


constructSet <- function(trioSet, CHR, id, states, ranges){
	open(baf(trioSet))
	open(logR(trioSet))
	i <- match(id[["O"]], offspringNames(trioSet))
	mads <- mad(trioSet)[i, ]
	S <- length(states)
	loglik <- array(NA, dim=c(2, nrow(trioSet), 3, S))
	dimnames(loglik) <- list(c("logR", "baf"),
				 featureNames(trioSet),
				 id,
				 states)
	object <- new("LikSet",
		      logR=as.matrix(logR(trioSet)[ ,i,]),
		      BAF=as.matrix(baf(trioSet)[ ,i , ]),
		      featureData=featureData(trioSet),
		      loglik=loglik)
	object$MAD <- mads
	fData(object)$range.index <- NA
	ir1 <- IRanges(start=position(object), end=position(object))
	ir2 <- IRanges(start(ranges), end(ranges))
	mm <- findOverlaps(ir1, ir2)
	subject.index <- subjectHits(mm)
	## there should be no query that is in more than 1 subject
	qhits <- queryHits(mm)
	shits <- subjectHits(mm)
	fData(object)$range.index <- NA
	fData(object)$range.index[qhits] <- shits
	close(baf(trioSet))
	close(logR(trioSet))
	return(object)
}

shrinkTo <- function(x, x.0, DF.PRIOR){
	DF <- ncol(x)-1
	DF <- Ns-1
	DF[DF < 1] <- 1
	x.0 <- apply(x, 2, median, na.rm=TRUE)
	x <- (x*DF + x.0*DF.PRIOR)/(DF.PRIOR + DF)
	for(j in 1:ncol(x)) x[is.na(x[, j]), j] <- x.0[j]
	return(x)
}


dna <- function(object) harmonizeDnaLabels(phenoData2(object[[1]])[, "DNA.Source", ])
plate <- function(object) phenoData2(object[[1]])[, "Sample.Plate", ]


fillInMissing <- function(rangeIndex){
	if(!any(is.na(rangeIndex))) return(rangeIndex)
	if(sum(is.na(rangeIndex)) > 1000 & length(unique(rangeIndex[!is.na(rangeIndex)])) == 1){
		## for calculating the posterior for a single range
		## essentially, ignoring what comes before and what comes after
		ii <- range(which(!is.na(rangeIndex)))
		if(ii[1] > 1){
			rangeIndex[1:(ii[1]-1)] <- 0
		}
		if(ii[2] < length(rangeIndex)){
			rangeIndex[(ii[2]+1):length(rangeIndex)] <- 2
		}
	} else {
		ii <- which(is.na(rangeIndex))
		if(max(ii) < length(rangeIndex)){
			rangeIndex[ii] <- rangeIndex[ii+1]
		} else{
			rangeIndex[max(ii)] <- rangeIndex[max(ii)-1]
			if(any(ii != max(ii))){
				iii <- ii[ii != max(ii)]
				rangeIndex[iii] <- rangeIndex[iii+1]
			}
		}
	}
	return(rangeIndex)
}



mosaicProb <- function(bf, sd.mosaic, sd0, sd.5, sd1, rangeIndex, normalCN=FALSE, object){
	initialP <- 0.99## initial probabilty not mosaic
	iter <- 1
	##ri=range.index(object)
	while(any(is.na(rangeIndex))){
		rangeIndex <- fillInMissing(rangeIndex)
		iter <- iter+1
		if(iter > 5) break()
	}
	f1 <- tnorm(bf, 0, sd0)
	f2 <- tnorm(bf, 1, sd1)
	f3 <- tnorm(bf, 0.25, sd.mosaic, lower=0.05, upper=0.4)
	f4 <- tnorm(bf, 0.75, sd.mosaic, lower=0.6, upper=0.95)
	if(normalCN) f5 <- tnorm(bf, 0.5, sd.5)
	tau <- matrix(initialP, nrow(f1), ncol(f1))
	LL <- sapply(split(rangeIndex, rangeIndex), length)
	p.out <- 1e-10
	## updating f3,f4 might help
	for(k in 1:3){
		if(!normalCN){
			T1.num <- tau * ((1-p.out)*(0.5*f1 + 0.5*f2) + p.out)
			T.den <- tau * ((1-p.out)*(0.5*f1 + 0.5*f2) + p.out) + (1-tau)*((1-p.out)*0.25*f1 + 0.25*f3 + 0.25*f4 + 0.25*f2) + p.out
		} else {
			T1.num <- tau * ((1-p.out)*(1/3*f1 + 1/3*f2 + 1/3*f5) + p.out)
			T.den <- tau * ((1-p.out)*(1/3*f1 + 1/3*f2 + 1/3*f5) + p.out) + (1-tau)*((1-p.out)*(0.25*f1 + 0.25*f3 + 0.25*f4 + 0.25*f2) + p.out)
		}
		T1 <- T1.num/T.den
		tau.F <- sapply(split(T1[, 1], rangeIndex), mean, na.rm=TRUE)
		tau.M <- sapply(split(T1[, 2], rangeIndex), mean, na.rm=TRUE)
		tau.O <- sapply(split(T1[, 3], rangeIndex), mean, na.rm=TRUE)
		tau.F <- rep(tau.F, LL)
		tau.M <- rep(tau.M, LL)
		tau.O <- rep(tau.O, LL)
		tau.next <- cbind(tau.F, tau.M, tau.O)
		if(abs(sum(tau.next - tau, na.rm=TRUE)) < 1) break()
		tau <- tau.next
		##tau.next <- apply(T1, 2, mean, na.rm=TRUE)
		##tau <- tau.next
	}
	return(tau)
}




emissionB <- function(p.roh=0.01, q.roh=1-p.roh,
		      p.mos=0.01, q.mos=1-p.mos,
		      p.out=1e-15, q.out=1-p.out,
		      sd0, sd1, sd.5,
		      object){
	## model emission as a mixture of normals (genotype is correct) and a uniform (error)
	## Wang et al. use mixture probabilities from a binomial
	##
	##   pi ~ binomial(C(z), allele freq)
	##   and integrates out the number of copies of the B allele.
	##
	##   Below, I,ve just used a mixture model.  I have not integrated out the
	##      copy number of the B allele, nor do I make use of MAF estimates.
	##   TN vs dnorm shouldn't make much diff.
	p.mosaic=p.mos
	sd.mosaic <- 0.05
	bf <- baf(object)
	if(FALSE){
		ii <- which(range.index(object)==i)
		q.mos <- mosaicProb(bf=bf, sd.mosaic=sd.mosaic/2,
				       sd=sd0, sd.5=sd.5, sd1=sd1,
				       rangeIndex=range.index(object),
				       normalCN=FALSE)
		p.mos <- 1-q.mos
	}
	##p.mos <- prMosaic
	q.mos <- 1-p.mos
	p.out <- 1e-10
	q.out <- p.out
	TN0 <- tnorm(bf, 0, sd0);
	TN1 <- tnorm(bf, 1, sd1)
	TN.3 <- tnorm(bf, 1/3, sd.5)
	TN.6 <- tnorm(bf, 2/3, sd.5)
	TN.25 <- tnorm(bf, 0.25, sd.5);
	TN.5 <- tnorm(bf, 0.5, sd.5);
	TN.75 <- tnorm(bf, 0.75, sd.5)
	TN.M1 <- tnorm(bf, 0.25, sd.mosaic, lower=0.02, upper=0.46);
	TN.M2 <- tnorm(bf, 0.75, sd.mosaic, lower=0.54, upper=0.98);
	## p1 <- dunif(bf, 0, 1)  if you do a uniform for mosaic, it allows baf's near 0.5 to creep in
	p1 <- 0.5*TN.M1 + 0.5*TN.M2
	beta.hemizygous <- (1-q.out)*(p.mos*p1 + q.mos*(0.5*TN0+0.5*TN1)) + p.out
	q.mos2 <- mosaicProb(bf=bf, sd.mosaic=sd.mosaic,
			    sd=sd0, sd.5=sd.5, sd1=sd1,
			    rangeIndex=range.index(object),
			    normalCN=TRUE)
	q.mos2[is.nan(q.mos2)] <- 0.99
	p.mos2 <- 1-q.mos2
	pr1 <- p1
	pr2 <- 0.5*TN0 + 0.5*TN1
	pr3 <- p.roh*(p.mos2*pr1 + q.mos2*pr2)
	pr4 <- p1
	pr5 <- 1/3*TN0 + 1/3*TN.5 + 1/3*TN1
	pr6 <- q.roh*(p.mos2*pr4 +  q.mos2*pr5)
	beta.normal <- pr3+pr6
	loglik(object)["baf", , , 1] <- dunif(bf, 0, 1)
	loglik(object)["baf", , , 2] <- beta.hemizygous
	loglik(object)["baf", , , 3] <- beta.normal
	loglik(object)["baf", , , 4] <- 1/4*TN0 + 1/4*TN.3 + 1/4*TN.6 + 1/4*TN1
	loglik(object)["baf", , , 5] <- 1/5*TN(bf, 0, sd0) + 1/5*TN(bf, 1/4, sd.5) + 1/5*TN(bf, 0.5, sd.5) + 1/5*TN(bf, 0.75, sd.5) + 1/5*TN(bf, 1, sd1)
	loglik(object)["baf", , , ] <- log(loglik(object)["baf", , , ])
	if(FALSE){
		LLB <- loglik(object)["baf", range.index(object)==i , , ]
		LLR <- loglik(object)["logR", range.index(object)==i , , ]
		ii <- which(range.index(object)==i)
		b <- baf(object)[ii, 3]
		ii2 <- which(b > 0.45 & b < 0.55)
		apply(LLB[, 1, ], 2, sum, na.rm=TRUE)
		apply(LLR[, 1, ], 2, sum, na.rm=TRUE)
		apply(LLB[, 2, ], 2, sum, na.rm=TRUE)
		apply(LLB[, 3, ], 2, sum, na.rm=TRUE)
		apply(LLR[, 3, ], 2, sum, na.rm=TRUE)
		LL <- LLR + LLB
		LLT <- matrix(NA, 3, 5)
		for(j in 1:3) LLT[j, ] <- apply(LL[, j, ], 2, sum, na.rm=TRUE)
		rownames(LLT) <- c("F", "M", "O")
		colnames(LLT) <- paste("CN_", 0:4, sep="")
	}
	return(object)
}

	rohProb <- function(x, sd0, sd.5, sd1, rangeIndex){
		initialP <- 0.01## initial probabilty for region of homozygosity
		fAA <- tnorm(x, 0, sd0)
		fBB <- tnorm(x, 1, sd1)
		fAB <- tnorm(x, 0.5, sd.5)
		tau <- matrix(initialP, nrow(fAA), ncol(fAA))
		LL <- sapply(split(rangeIndex, rangeIndex), length)
		p.out <- 1e-10
		## updating f3,f4 might help
		for(k in 1:3){
			T1.num <- tau * ((1-p.out)*(0.5*fAA + 0.5*fBB) + p.out)
			T.den <- tau * ((1-p.out)*(0.5*fAA + 0.5*fBB) + p.out) + (1-tau)*((1-p.out)*(1/3*fAA + 1/3*AB + 1/3*fBB) + p.out)
			T1 <- T1.num/T.den
			tau.F <- sapply(split(T1[, 1], rangeIndex), mean, na.rm=TRUE)
			tau.M <- sapply(split(T1[, 2], rangeIndex), mean, na.rm=TRUE)
			tau.O <- sapply(split(T1[, 3], rangeIndex), mean, na.rm=TRUE)
			tau.F <- rep(tau.F, LL)
			tau.M <- rep(tau.M, LL)
			tau.O <- rep(tau.O, LL)
			tau.next <- cbind(tau.F, tau.M, tau.O)
			if(abs(sum(tau.next - tau, na.rm=TRUE)) < 1) break()
			tau <- tau.next
		}
		tau
	}

emissionLR <- function(mu.logr, CN.MIN, CN.MAX, prMosaic,
		       prOutlier, sds, object, ranges){
	##CN.MAX=10
	##CN.MIN=-10
	p.out <- prOutlier
	q.out <- 1-p.out
	lR <- logR(object)
	if(any(lR < CN.MIN, na.rm=TRUE)) lR[lR < CN.MIN] <- CN.MIN
	if(any(lR > CN.MAX, na.rm=TRUE)) lR[lR > CN.MAX] <- CN.MAX
	UNIF <- dunif(lR, CN.MIN, CN.MAX)
	##pRoh <- rohProb(x=baf(object), sd0=sd0, )
##	bf <- baf(object)
##	ri <- range.index(object)
##	LL <- sapply(split(ri,ri), length)
##	phom <- function(x) mean(x < 0.02 | x > 0.98, na.rm=TRUE)
##	pHet.F <- sapply(split(bf[,1], ri), phom)
##	pHet.F <- rep(pHet.F, LL)
##	pHet.M <- sapply(split(bf[,2], ri), phom)
##	pHet.M <- rep(pHet.M, LL)
##	pHet.O <- sapply(split(bf[,3], ri), phom)
##	pHet.O <- rep(pHet.O, LL)
##	mu <- matrix(mu.logr[2], nrow(lR), ncol(lR))
##	mu[, 1] <- ifelse(pHet.F < 0.01, mu[,1]/3,  mu[,1])
##	mu[, 2] <- ifelse(pHet.M < 0.01, mu[,1]/3,  mu[,2])
##	mu[, 3] <- ifelse(pHet.O < 0.01, mu[,1]/3,  mu[,3])
	##MU <- matrix(-0.5, nrow(lR), ncol(lR))
	##MU2 <- MU * (1-p
	loglik(object)["logR", , , 1] <- q.out * dunif(lR, CN.MIN, -1) + p.out*UNIF
	##propGreaterZero <- sapply(split(lR, range.index(object)), function(x) mean(x>0,na.rm=T))
	##test whether = 0.5
	##p <- ifelse(propGreaterZero > 0.4, prMosaic=0.01, 0.99)
 	##prMosaic <- ifelse(p > 0.4, 0.01, 1-
	## estimate prMosaic
	## perhaps EM with 2 or 3 iterations to keep it fast
##	loglik(object)["logR", , , 2] <- q.out * (prMosaic * dnorm(lR, mu.logr[2]/2, sds) + (1-prMosaic) * dnorm(lR, mu.logr[2], sds)) + p.out*UNIF
	loglik(object)["logR", , , 2] <- q.out * (prMosaic * dnorm(lR, mu.logr[2]/2, sds) + (1-prMosaic) * dnorm(lR, mu.logr[2], sds)) + p.out*UNIF
	loglik(object)["logR", , , 3] <- q.out * dnorm(lR, mu.logr[3], sds) + p.out*UNIF
	loglik(object)["logR", , , 4] <- q.out * dnorm(lR, mu.logr[4], sds) + p.out*UNIF
	loglik(object)["logR", , , 5] <- q.out * dnorm(lR, mu.logr[5], sds) + p.out*UNIF
	loglik(object)["logR", , , ] <- log(loglik(object)["logR", , , ])
	if(FALSE){
		q.out[ii, ]
		LLB <- loglik(object)["baf", range.index(object)==i , , ]
		LLR <- loglik(object)["logR", range.index(object)==i , , ]
		ii <- which(range.index(object)==i)
		apply(LLR[, 1, ], 2, sum, na.rm=TRUE)
		apply(LLR[, 2, ], 2, sum, na.rm=TRUE)
		apply(LLR[, 3, ], 2, sum, na.rm=TRUE)
		LL <- LLR + LLB
		LLT <- matrix(NA, 3, 5)
		for(j in 1:3) LLT[j, ] <- apply(LL[, j, ], 2, sum, na.rm=TRUE)
		rownames(LLT) <- c("F", "M", "O")
		colnames(LLT) <- paste("CN_", 0:4, sep="")
	}
	return(object)
}

computeLoglik <- function(id,
			  trioSet, #L,
			  ranges,
			  mu.logr=c(-2, -0.5, 0, 0.3, 0.75),
			  states=0:4,
			  baf.sds=c(0.02, 0.03, 0.02),
			  ##baf.sds=c(0.015, 0.02, 0.015),
			  ##prGtCorrect=0.999, ##prob genotype is correct
			  prOutlier.logR=0.01,
			  prOutlier.baf=1e-5,
			  prMosaic=0.01,
			  df0=50,
			  returnEmission=FALSE){
	CHR <- chromosome(trioSet)[1]
	##p1 <- prGtCorrect; rm(prGtCorrect)
	p1 <- 1-prOutlier.logR
	## one obvious thing that p1 could depend on is the
	## minor allele frequency.  If rare, p1 is smaller
	stopifnot(all(!is.na(match(id, s(fullId(trioSet))))))
	object <- constructSet(trioSet, CHR, id, states=states, ranges=ranges)
	j <- match(id[["O"]], offspringNames(trioSet))
	sds.sample <- mad(trioSet)[j, ]
	stopifnot(all(!is.na(sds.sample)))
	sds.sample <- matrix(sds.sample, nrow(trioSet), 3, byrow=TRUE)
	open(logR(trioSet))
##	sds.marker <- rowMAD(logR(trioSet)[, , "O"], na.rm=TRUE)
	sds.marker <- fData(trioSet)$marker.mad  ## these can be really big in CNV
	##1. shrink the marker to the marker sds
	##   - at CNV the marker estimates will be much too high
	##   - could be too small for markers that don't really work
	df1 <- nrow(phenoData(trioSet))/3
	df2 <- min(length(sds.marker)-1, 10e3)
	sds.marker <- (df1*sds.marker + df2*median(sds.marker,na.rm=TRUE))/(df1+df2)
	sds.marker <- matrix(sds.marker, nrow(object), 3, byrow=FALSE)
	##2. shrink to the sample-level variance
	## -- how much do we want to shrink towards a sample-level estimate of noise?
	sds <- (sds.marker * df1 + (df1*4)*sds.sample)/(df1 + df1*4)
	lR <- logR(object)

	bf <- baf(object)
	if(!is(baf.sds, "array")){
		stopifnot(length(baf.sds)==3)
		sd0 <- baf.sds[1]
		sd.5 <- baf.sds[2]
		sd1 <- baf.sds[3]
	} else{
		sample.index <- match(id[3], rownames(baf.sds))
		stopifnot(length(sample.index)==1)
		baf.sds2 <- baf.sds[sample.index, , ]
		nr <- nrow(bf)
		sd0 <- matrix(baf.sds2[, "AA"], nr, 3, byrow=TRUE)
		sd.5 <- matrix(baf.sds2[, "AB"], nr, 3, byrow=TRUE)
		sd1 <- matrix(baf.sds2[, "BB"], nr, 3, byrow=TRUE)
	}
	q.mosaic <- mosaicProb(bf=bf,
			       sd.mosaic=0.05,
			       sd0=sd0, sd.5=sd.5,
			       sd1=sd1,
			       rangeIndex=range.index(object),
			       object=object)
	q.mosaic[is.nan(q.mosaic)] <- 0.99 ## prob not mosaic
	p.mosaic <- 1-q.mosaic
	## the uniform needs to cover the support
	CN.MIN <- -5; CN.MAX <- 1.5
	object <- emissionLR(mu.logr=mu.logr, CN.MIN=CN.MIN, CN.MAX=CN.MAX,
			     prMosaic=p.mosaic,
			     prOutlier=prOutlier.logR,
			     sds=sds,
			     object=object,
			     ranges=ranges)
	##p1 <- 1-prOutlier.baf

	object <- emissionB(p.roh=0.01,
			    p.mos=p.mosaic,
			    p.out=prOutlier.baf,
			    sd0=sd0, sd1=sd1, sd.5=sd.5,
			    object=object)
##	p.mos <- prMosaic
##	q.mos <- 1-p.mos
##	p.out <- 1e-15
##	q.out <- p.out
##	TN0 <- TN(bf, 0, sd0);
##	TN1 <- TN(bf, 1, sd1)
##	TN.3 <- dnorm(bf, 1/3, sd.5)
##	TN.6 <- dnorm(bf, 2/3, sd.5)
##	TN.25 <- dnorm(bf, 0.25, sd.5);
##	TN.5 <- dnorm(bf, 0.5, sd.5);
##	TN.75 <- dnorm(bf, 0.75, sd.5)
##	loglik(object)["baf", , , 1] <- dunif(b, 0, 1)
##	loglik(object)["baf", , , 2] <- p.mos*dunif(b, 0, 1) + q.mos*(0.5*TN0 + 0.5*TN1)
##	loglik(object)["baf", , , 3] <- q.roh*(p.mos*dunif(b, 0, 1) + q.mos*(1/3*TN0 + 1/3*TN.5 + 1/3*TN1)) +
##		p.roh*(p.mos*dunif(b, 0, 1) + q.mos *(0.5*TN0 + 0.5*TN1))
##	loglik(object)["baf", , , 4] <- p.mos*dunif(b, 0, 1) + q.mos*(1/4*TN0 + 1/4*TN.3 + 1/4*TN.6 + 1/4*TN1)
##	loglik(object)["baf", , , 5] <- p.mos*dunif(b, 0, 1) + q.mos*(1/5*TN(bf, 0, sd0) + 1/5*TN(bf, 1/4, sd.5) + 1/5*TN(bf, 0.5, sd.5) + 1/5*TN(bf, 0.75, sd.5) + 1/5*TN(bf, 1, sd1))
##	return(object)
##}
##
##	ploh <- 0.01
##	p.m <- prMosaic
##	q.m <- 1-p.m
##	pi.roh <- 0.01 ## region of homozygosity
##	loglik(object)["baf", , , 1] <-  1
##	##loglik(object)["baf", , , 2] <- p1*((1/2*TN(bf, 0, sd0) + 1/2*TN(bf, 1, sd1))) + (1-p1)  ## * dunif(bf, 0, 1) = 1
##	t0 <- TN(bf, 0, sd0); t1 <- TN(bf, 1, sd1)
##	t0.25 <- dnorm(bf, 0.25, sd.5); t0.75 <- dnorm(bf, 0.75, sd.5)
####	tmp <- p1*(
####					    p2*(1/2*t0 + 1/2*t1) + (1-p2)*(1/4*dunif(bf, 0, 0.2) + 1/4*dunif(bf, 0.8,1) + 1/4*t0.25 + 1/4*t0.75)
####					    ) + (1-p1)  ## * dunif(bf, 0, 1) = 1
##
##	tmp <- q.m*(0.5*t0 + 0.5*t1) + p.m*dunif(bf, 0, 1)
##					    p2*(1/2*t0 + 1/2*t1) + (1-p2)*(1/2*dunif(bf, 0, 0.4) + 1/2*dunif(bf, 0.6,1))
##					    ) + (1-p1)  ## * dunif(bf, 0, 1) = 1
####	tmp <- p1*(
####					    p2*(1/2*t0 + 1/2*t1) + (1-p2)*(1/2*dunif(bf, 0, 0.4) + 1/2*dunif(bf, 0.6,1))
####					    ) + (1-p1)  ## * dunif(bf, 0, 1) = 1
##	loglik(object)["baf", , , 2] <- tmp
##	## a better solution for loh would be to add a latent indicator for
##	## LOH, and then integrate this out of the likelihood
##
##	tmp1 <- p1*(p2*(1/3*TN(bf, 0, sd0) + 1/3*TN(bf, 0.5, sd.5) + 1/3*TN(bf, 1, sd1)) + (1-p2)*(1/3*dunif(bf, 0, 0.2) + 1/3*dunif(bf, 0.3, 0.7) + 1/3*dunif(bf, 0.8, 1)))+ (1-p1)
##	##tmp2 <- p1*(p2*(1/2*t0 + 1/2*t1) + (1-p2)*(1/2*dunif(bf, 0, 0.2) + 1/2*dunif(bf,0.8,1))) + 1-p1
##	Lik.Nor <- matrix(NA, nrow(tmp1), ncol(tmp1))
##	ri <- range.index(object)
##	if(mean(is.na(ri)) > 0.50){
##		##
##		## For computing the bayes factor for just 1 range (e.g., a
##		## range identified by another method)
##		##
##		##   -- would be better to explicitly indicate this
##		##
##		## assign 'fake' range indices to the other markers
##		first <- which(!is.na(ri))[1]
##		last <- rev(which(!is.na(ri)))[1]
##		stopifnot(all(!is.na(ri[first:last])))
##		ri[1:(first-1)] <- ri[first]-1
##		ri[(last+1):length(ri)] <- ri[last]+1
##	}
##	counter <- 1
##	while(any(is.na(ri))){
##		isna.index <- which(is.na(ri))
##		if(any(isna.index == length(ri))){
##			ri[length(ri)] <- ri[length(ri)-1]
##		}
##		isna.index <- which(is.na(ri))
##		if(length(isna.index) > 0)
##			ri[isna.index] <- ri[isna.index+1]
##		counter <- counter+1
##		if(counter > 5) stop("problems with nas in range indx")
##	}
##	for(l in 1:3){
##		lik.normal <- sapply(split(log(tmp1[, l]), ri), sum,na.rm=T)
##		lik.loh <- sapply(split(log(tmp[, l]), ri), sum, na.rm=T)
##		isloh <- lik.loh > lik.normal
##		f <- sapply(split(ri, ri), length)
##		isloh <- rep(isloh,f)
##		Lik.Nor[, l] <- tmp1[, l]*(1-isloh) + isloh*(tmp[, l])
##	}
##	##loglik(object)["baf", , , 3] <- p1*((1-ploh)*(1/3*TN(bf, 0, sd0) + 1/3*TN(bf, 0.5, sd.5) + 1/3*TN(bf, 1, sd1)) + ploh*(1/2*t0 + 1/2*t1))+ (1-p1)
##	loglik(object)["baf", , , 3] <- Lik.Nor
##	loglik(object)["baf", , , 4] <- p1*((1/4*TN(bf, 0, sd0) + 1/4*TN(bf, 1/3, sd.5) + 1/4*TN(bf, 2/3, sd.5) + 1/4*TN(bf, 1, sd1))) + (1-p1)
##	loglik(object)["baf", , , 5] <- p1*((1/5*TN(bf, 0, sd0) + 1/5*TN(bf, 1/4, sd.5) + 1/5*TN(bf, 0.5, sd.5) + 1/5*TN(bf, 3/4, sd.5) + 1/5*TN(bf, 1, sd1))) + (1-p1)
##	loglik(object)["baf", , , ] <- log(loglik(object)["baf", , , ])
	if(returnEmission){
		return(object)
	}
	if(FALSE){
		LLR <- loglik(object)["logR", range.index(object)==i, ,  ]
		LLB <- loglik(object)["baf", range.index(object)==i , , ]
		ii <- which(range.index(object)==i)
		b=bf[ii, 3]
		##colnames(tmp4) <- c("hemizygous", "normal", "log r")
		##loglik(object)["logR", , , 3] <- p1 * dnorm(lR, mu.logr[3], sds) + (1-p1)*UNIF
		r <- cbind(LLR[, 3, 2:3], lR[ii, 3])
		##p1 <- prOutlier[1]
		##tmp <- p1 * dnorm(lR[ii, 1], mu.logr[3], sds[ii, 1]) + (1-p1)*UNIF[ii,1]
		apply(LLB[, 1, ], 2, sum, na.rm=TRUE)
		apply(LLR[, 1, ], 2, sum, na.rm=TRUE)
		apply(LLB[, 2, ], 2, sum, na.rm=TRUE)
		apply(LLR[, 2, ], 2, sum, na.rm=TRUE)
		apply(LLB[, 3, ], 2, sum, na.rm=TRUE)
		apply(LLR[, 3, ], 2, sum, na.rm=TRUE)
		LLT <- matrix(NA, 3, 5)
		LL <- LLR + LLB
		LLT <- matrix(NA, 3, 5)
		for(j in 1:3) LLT[j, ] <- apply(LL[, j, ], 2, sum, na.rm=TRUE)
		rownames(LLT) <- c("F", "M", "O")
		colnames(LLT) <- paste("CN_", states, sep="")
	}
	if(FALSE){
		pdf("../../figures/mixture.pdf", width=8, height=6)
		par(las=1, mfrow=c(1,1))
		plot(CN.MIN:CN.MAX,
		     type="n", xlab="", ylab="", ylim=c(0,3),
		     yaxs="i", xlim=c(-2, 1.5))
		##
		##
		colors <- c("grey0",
			    "grey20",
			    "grey40",
			    "grey80",
			    "white", "lightblue")
		colors <- makeTransparent(colors, alpha=0.3)
		## homozygous deletion
		rect(xleft=-4, xright=-1,
		     ybottom=0, ytop=1/(-1-CN.MIN),
		     col=colors[1])
		##
		## hemizygous deletion
		p <- seq(0, 1, by=0.001)
		x <- qnorm(p, mean=mu.logr[2], sds[1,1])
		##plot(function(x) dnorm(x, log=TRUE), -60, 50,
		##     main = "log { Normal density }")
		y <- dnorm(x, mean=mu.logr[2], sd=sds[1,1])
		##y <- y/max(y[is.finite(y)])
		polygon(c(x, rev(x)), c(rep(0, length(x)), rev(y)), col=colors[2])
		rect(xleft=-1, xright=0,
		     ybottom=0, ytop=1,
		     col=colors[2])
		## normal
		x <- qnorm(p, mean=mu.logr[3], sds[1,1])
		y <- dnorm(x, mean=mu.logr[3], sd=sds[1,1])
		polygon(c(x, rev(x)), c(rep(0, length(x)), rev(y)), col=colors[3])
		## amplification
		x <- qnorm(p, mean=mu.logr[4], sds[1,1])
		y <- dnorm(x, mean=mu.logr[4], sd=sds[1,1])
		polygon(c(x, rev(x)), c(rep(0, length(x)), rev(y)), col=colors[4])
		## 2 copy dup
		x <- qnorm(p, mean=mu.logr[5], sds[1,1])
		y <- dnorm(x, mean=mu.logr[5], sd=sds[1,1])
		polygon(c(x, rev(x)), c(rep(0, length(x)), rev(y)), col=colors[5])
		## outlier
		rect(xleft=CN.MIN, xright=CN.MAX,
		     ybottom=0, ytop=1/(CN.MAX-CN.MIN), col=colors[6])
		legend("left", bty="n", fill=colors,
		       legend=c("homozygous del",
		       "hemizygous del",
		       "normal",
		       "1 copy dup",
		       "2 copy dup", "outlier"))
		dev.off()

		par(las=1, mfrow=c(1,1))
		plot(0:1,
		     type="n", xlab="", ylab="", ylim=c(0,80),
		     yaxs="i", xlim=c(0, 1))
		##
		##
		colors <- c("grey0",
			    "grey20",
			    "grey40",
			    "grey80",
			    "white", "lightblue")
		colors <- makeTransparent(colors, alpha=0.3)
		## homozygous deletion
		rect(xleft=0, xright=1,
		     ybottom=0, ytop=1,
		     col=colors[6])
		## hemizygous deletion
		require(msm)
		p <- seq(0, 1, by=0.001)
		x <- qtnorm(p, mean=0, 0.02, lower=0, upper=1)
		##plot(function(x) dnorm(x, log=TRUE), -60, 50,
		##     main = "log { Normal density }")
		y <- dtnorm(x, mean=0, 0.02, lower=0, upper=1)
		##y <- y/max(y[is.finite(y)])
		x[1]=0
		polygon(c(x, rev(x)), c(rep(0, length(x)), rev(y)), col=colors[2])
		x <- qtnorm(p, mean=1, 0.02, lower=0, upper=1)
		y <- dtnorm(x, mean=1, 0.02, lower=0, upper=1)
		##y <- y/max(y[is.finite(y)])
		x[length(x)]=1
		x[1]=0
		polygon(c(x, rev(x)), c(rep(0, length(x)), rev(y)), col=colors[2])
		## allowing for mosaic cn
		rect(0, 0, 0.5, 0.5, col=colors[2])
		rect(0.5, 0, 1, 0.5, col=colors[2])
		## normal
		x <- qnorm(p, mean=0.5, 0.03)
		y <- dnorm(x, mean=0.5, 0.03)
		##y <- y/max(y[is.finite(y)])
		x[length(x)]=1
		x[1]=0
		polygon(c(x, rev(x)), c(rep(0, length(x)), rev(y)), col=colors[3])
		## single dup
		x <- qnorm(p, mean=1/3, 0.03)
		y <- dnorm(x, mean=1/3, 0.03)
		##y <- y/max(y[is.finite(y)])
		x[length(x)]=1
		x[1]=0
		polygon(c(x, rev(x)), c(rep(0, length(x)), rev(y)), col=colors[4])
		x <- qnorm(p, mean=2/3, 0.03)
		y <- dnorm(x, mean=2/3, 0.03)
		##y <- y/max(y[is.finite(y)])
		x[length(x)]=1
		x[1]=0
		polygon(c(x, rev(x)), c(rep(0, length(x)), rev(y)), col=colors[4])
		## 2 copy duplication
		##...
		legend("left", bty="n", fill=colors,
		       legend=c("homozygous del",
		       "hemizygous del",
		       "normal",
		       "1 copy dup",
		       "2 copy dup", "outlier"))



	}
	return(object)
}

rowMAD <- function(x, y, ...){
	##notna <- !is.na(x)
	mad <- 1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)
	return(mad)
}

trioStates <- function(states=0:4){
	trio.states <- as.matrix(expand.grid(states, states, states))
	index <- which(trio.states[, 1] == 0 & trio.states[, 2] == 0 & trio.states[, 3] > 0)
	trio.states <- trio.states+1
	colnames(trio.states) <- c("F", "M", "O")
	## 125 possible
	## remove 00 > 0 as possibilities
	trio.states <- trio.states[-index, ]
}

trioStateNames <- function(trio.states){
	if(missing(trio.states)) trio.states <- trioStates(0:4)
	paste(paste(trio.states[,1], trio.states[,2], sep=""), trio.states[,3], sep="")
}


transitionProbability <- function(states=0:4, epsilon=1-0.999){
	off.diag <- epsilon/(length(states)-1)
	tpm <- matrix(off.diag, length(states), length(states))
	diag(tpm) <- 1-epsilon
	tpm
}

initialStateProbs <- function(states=0:4, normal.index=3, epsilon=0.01){
	initial.state.probs <- rep(epsilon/(length(states)-1), length(states))
	initial.state.probs[normal.index] <- 1-epsilon
	initial.state.probs
}

readTable1 <- function(states=0:4, a=0.0009){
	S <- length(states)
	tmp <- array(NA, dim=rep(S,3))
	dimnames(tmp) <- list(paste("F", states, sep=""),
			      paste("M", states, sep=""),
			      paste("O", states, sep=""))
	tmp["F0", "M0", ] <- c(1, rep(0,4))
	tmp["F0", "M1", ] <- c(0.5, 0.5, 0, 0, 0)
	tmp["F0", "M2", ] <- c(0.5*a, 1-a, 0.5*a, 0, 0)
	tmp["F0", "M3", ] <- c(0.5*a, 0.5*(1-a), 0.5*(1-a), 0.5*a, 0)
	tmp["F0", "M4", ] <- c(0, 0.25, 0.5, 0.25, 0)

	tmp["F1", "M1", ] <- c(0.25, 0.5, 0.25, 0, 0)
	tmp["F1", "M2", ] <- c(0.25, 0.5-0.25*a, 0.5-0.25*a, 0.25*a, 0)
	tmp["F1", "M3", ] <- c(0.25*a, 0.25, 0.5*(1-a), 0.25, 0.25*a)
	tmp["F1", "M4", ] <- c(0, 0.125, 0.375, 0.375, 0.125)

	tmp["F2", "M2", ] <- c(0.25*a^2, a*(1-a), (1-a)^2 + 0.5*a^2, a*(1-a), 0.25*a^2)
	tmp["F2", "M3", ] <- c(0.25*a^2, 0.75*a*(1-a), 0.5*(1-a)^2+0.25*a*(1-a)+0.25*a^2,
			       0.5*(1-a)^2+0.25*a*(1-a)+0.25*a^2,
			       0.75*a*(1-a)+0.25*a^2)
	tmp["F2", "M4", ] <- c(0,0.125*a, 0.25, 0.5-0.25*a,0.25+0.125*a)

	tmp["F3", "M3", ] <- c(0.25^2, 0.5*a*(1-a), 0.5*a*(1-a)+0.25*(1-a)^2, 0.5*(1-a)^2+0.5^2,
			       0.25*(1-a)^2 + a*(1-a)+0.25^2)
	tmp["F3", "M4", ] <- c(0, 0.125*a, 0.125*(1+a), 0.125*(3-2*a), 0.5)

	tmp["F4", "M4", ] <- c(0, 0, 0.0625, 0.25, 0.6875)
	return(tmp)
}

lookUpTable1 <- function(table1, state){
	if(is.na(table1[state[1], state[2], state[3]])){
		return(table1[state[2], state[1], state[3]])
	} else {
		return(table1[state[1], state[2], state[3]])
	}
}

readTable3 <- function(a=0.009){
	## initialize with small value to avoid -Inf
	results <- .C("calculateCHIT", a=a, M=array(0, dim=c(rep(5,6))))$M
	## Make sure to transpose!
	aperm(results)
}

lookUpTable3 <- function(table3, state.prev, state.curr){
	f1 <- state.prev[1]
	f2 <- state.curr[1]
	m1 <- state.prev[2]
	m2 <- state.curr[2]
	o1 <- state.prev[3]
	o2 <- state.curr[3]
	return(table3[f1, f2, m1, m2, o1, o2])
}

## p.00 = Prob(S_l, S_l-1 | I_l-1 = 0, I_l = 0)
##      = lookuptable3        ## both Mendelian
## p.10 = Prob(S_l, S_l-1 | I_l-1 = 1, I_l = 0) ## previous is nonmendelian
##      = 1/5*Prob(S_l | s_l,m, s_l,f, I_l=0)
##      = 1/5*lookUpTable(current state)
## p.01 = Prob(S_l, S_l-1 | I_l-1 = 0, I_l = 1)
##      = Prob(S_l-1 | s_l-1,m, s_l-1,f, ...)
##      = 1/5*lookUpTable(previous state)  ##previous is mendelian
## p.11 = Prob(S_l, S_l-1 | I_l-1 = 1, I_l = 1)
##      = Prob(s_l | s_l-1)* Prob(s_l-1 | I_l-1=1)
##      = tau.o * 1/5, where
##  tau.o = Prob(s_l | s_l-1 )
## tauI.00 = Prob(I_l = 0 | I_l-1 = 0)
## tauI.10 = Prob(I_l = 0 | I_l-1 = 1)
## tauI.01 = Prob(I_l = 1 | I_l-1 = 0)
## tauI.11 = Prob(I_l = 1 | I_l-1 = 1)
## piI0 = Prob(I_l-1 = 0)
## piI1 = Prob(I_l-1 = 1)
## beta.r = Prob(r | ...) = P(r_f|...)*P(r_m|...)*P(r_o|...)
## beta.b = Prob(b | ...)
## beta = beta.r*beta.b
## tau.f  = Prob(s_l,f | s_l-1,f)
## tau.m  = Prob(s_l,m | s_l-1,m)
## pi.f = Prob(s_l-1,f | theta)
## pi.m = Prob(s_l-1,m | theta)
## then
jointProb <- function(segment.index, ## so that we can insert a browser for a specific segment
		      state,
		      state.prev,
		      prob.nonMendelian,
		      log.pi,
		      tau,
		      table1,
		      table3,
		      log.lik){
	pi.f <- exp(log.pi[state[1]])
	pi.m <- exp(log.pi[state[2]])
	tau.o <- tau[state.prev[3], state[3]]
	tau.m <- tau[state.prev[2], state[2]]
	tau.f <- tau[state.prev[1], state[1]]
	p.00 <- lookUpTable3(table3, state.prev, state.curr=state) ## both Mendelian
	p.10 <- 1/5*lookUpTable1(table1, state) ## previous non-Mendelian * current Mendelian
	p.01 <- lookUpTable1(table1, state.prev) * 1/5 ## previous Mendlian * current non-mendelian
	p.11 <- 1/5*tau.o
	piI0 <- 1-prob.nonMendelian
	p.NM <- piI1 <- prob.nonMendelian
	## setting this to a small value will favor '2,2,1' versus '3,3,1' (for example)
	## setting prob.nonMendelian smaller would not have an effect
	tauI.11 <- tauI.00 <- 1-.01
	tauI.10 <- tauI.01 <- 1-tauI.11
	##beta <- exp(log.lik)
	##pr.off <- p.00*tauI.00*piI0 + p.10*tauI.10*piI1 + p.01*tauI.01*piI0 +p.11*tauI.11*piI1
	pr.off <- p.NM*(p.11*tauI.11 + p.10*tauI.10) + (1-p.NM)*(p.00*tauI.00+p.01*tauI.01)
	log.lik <- sum(log.lik) + log(pr.off*tau.m*pi.m*tau.f*pi.f)
	## lik.f*lik.m*lik.o*pr.off*tau.m*pi.m*tau.f*pi.f
	##log.lik[3] <- log.lik[3]+log(pr.off)
	##log.lik[2] <- log.lik[2]+log(tau.m*pi.m)
	##log.lik[1] <- log.lik[1]+log(tau.f*pi.f)
}

joint1 <- function(LLT, ##object,
		   trio.states,
		   tau,
		   log.pi,
		   normal.index,
		   segment.index,
		   state.index,
		   table1,
		   table3,
		   is.denovo=FALSE,
		   prob.nonMendelian=1.5e-6,
		   denovo.prev=FALSE,
		   state.prev) {
	Prob.DN <- prob.nonMendelian
	state <- trio.states[state.index, ]
	fmo <- c(LLT[1, state[1]], LLT[2, state[2]], LLT[3, state[3]])
	if(segment.index == 1 | is.null(state.prev)){
		## assume Pr(z_1,f | lambda) = Pr(z_2,m | lambda) = pi
		## For offspring, we have Pr(z_1,o | z_1,f, z_1,m, DN=0, 1)
		##    or 1/5 if DN=1
		##
		## if DN is 0 (not devovo), then many of the hidden
		##  states should have essentially an epsilon
		##  probability of occurring.
		##log.Prob.DN <- ifelse(is.denovo, log(Prob.DN), log(1-Prob.DN))
		pi.offspring <- c(lookUpTable1(table1, state),  1/5)
		lpr.offspring <- log(pi.offspring[1] * (1-Prob.DN) + pi.offspring[2]+Prob.DN)
		##pi.offspring <- pi.offspring[[is.denovo+1]]
		log.pi <- c(log.pi[state[1]], ## father
			    log.pi[state[2]], ## mother
			    lpr.offspring)
			    ##log(pi.offspring))## offspring
		##fmo <- apply(fmo, 2, sum, na.rm=TRUE)
		fmo <- fmo + log.pi
		log.emit <- fmo
##		for(j in 1:2) fmo[j] <- fmo[j]+log.pi[state.index]
##		f <- sum(fmo[[1]])+log.pi[state.index]
##		m <- sum(fmo[[2]])+log.pi[state.index]
##		o <- sum(fmo[[3]])+log(pi.offspring)
		## p.E138 top left: Pr(DN_1 = 1 | model)
		##  -- the probability the first marker is in a denovo region
		##log.Prob.DN <- ifelse(is.denovo, log(Prob.DN), log(1-Prob.DN))
		##pr.offspring <-
		##fmo[3] <- fmo[3]+log.Prob.DN
	} else{
		log.emit <- jointProb(segment.index=segment.index,
				      state=state,
				  state.prev=state.prev,
				  prob.nonMendelian=prob.nonMendelian,
				  log.pi=log.pi,
				  tau=tau,
				  table1=table1,
				  table3=table3,
				  log.lik=fmo)
##		##** Note that when state is the 'normal.state'**
##		##    (state.index == normal.index)
##		##    tau is the probability of staying in the same state
##		##
##		##for k = normal.index, it would do the right thing
##		## prob. leaving normal state to state k
##		##fmo <- apply(fmo, 2, sum, na.rm=TRUE)
##		for(j in 1:2) fmo[j] <- fmo[j]+log(tau[state.prev[j], state[j]])
##		##f <- log(tau[state.prev[1], state[1]]) + sum(fmo[[1]])
##		##m <- log(tau[state.prev[2], state[1]]) + sum(fmo[[2]])
##		##if(denovo.prev){
##		##
##		## Previous non-mendelian
##		##
##		## current Mendelian
##		tabled.value <- lookUpTable1(table1, state)
##		p1 <- 1/5*tabled.value
##		## current non-Mendelian
##		p2 <- 1/5*tau[state.prev[3], state[3]]  ## eq. 9
##		##fmo[3] <- log(p1 + p2)
##		##fmo[3] <- (1/5*tau[state.prev[3], state[3]]*
##		##fmo[3] <- log(1/5) + log(tau[state.prev[3], state[3]]) + fmo[3]
##		##
##		## Previous Mendelian
##		##
##		## current mendelian
##		tabled.value <- lookUpTable3(table3, state.prev, state.curr=state)
##		p3 <- tabled.value
##		## current non-mendelian
##		p4 <- 1/5*lookUpTable1(table1, state) ## eq. 10
##		##
##		fmo[3] <- fmo[3]+log(p1+p2+p3+p4)
		##fmo[3] <- log(p1+p2)
		##fmo[3] <- log(tabled.value) + fmo[3]
		##if(denovo.prev & is.denovo){
			## Equation 9: Wang et al.
			## if previous marker was in a range that was called denovo
			## Pr(z_j,o, z_j-1,o | DN_j=1, DN_j-1=1)
			##  = Pr(z_j,o | z_j-1,o, DN_j=1, DN_j-1=1) * Pr(z_j-1,o|DN_j-1=1)
			##    2nd term is 1/5
			##    1st term is
			## Pr(CN state of offspring in previous segment | previous segment was denovo) = 1/5
##			fmo[3] <- log(1/5) + log(tau[state.prev[3], state[3]]) + fmo[3]
##		}
##		if(denovo.prev & !is.denovo){
##			## Equation 10: Wang et al.
##			## previous marker was in a range that was not called denovo
##			tabled.value <- lookUpTable1(table1, state)
##			fmo[3] <- log(1/5) + log(tabled.value) + fmo[3]
##		}
##		if(!denovo.prev & is.denovo){
##			## Equation 10
##			tabled.value <- lookUpTable1(table1, state.prev)
##			fmo[3] <- log(1/5) + log(tabled.value) + fmo[3]
##		}
##		if(!denovo.prev & !is.denovo){
##			## 25 x 25 x 25 table available in source code of software
##			tabled.value <- lookUpTable3(table3, state.prev, state.curr=state)
##			fmo[3] <- log(tabled.value) + fmo[3]
##		}
	}
	## RS: comment 4/29/2011
	## add the probability of transitioning back to the normal state
	##for(j in 1:3) fmo[j] <- fmo[j] + log(tau[state[j], 3])
	##f <- f+log(tau[state[1], 3])
	##m <- m+log(tau[state[2], 3])
	##o <- o+log(tau[state[3], 3])
	##res <- sum(fmo)
	res <- sum(log.emit)
	stopifnot(!is.na(res))
	##
	## we all need a transition proability for not denovo -> denovo -> not denovo
	## (Equation 11)
	##f.m.o <- sum(fmo, na.rm=TRUE)
	## prob. back to  normal state
	##f.m.o <- f.m.o+tau[state.index, normal.index]
	return(res)
}

joint4 <- function(trioSet,
		   ranges, ## all the ranges from one subject , one chromosome
		   states,
		   baf.sds,
		   THR=-50,
		   mu.logr=c(-2,-0.5, 0, 0.3, 0.75),
		   log.pi,
		   tau,
		   normal.index,
		   a=0.0009,
		   verbose=TRUE,
		   prOutlier=c(0.01, 1e-5),
		   prMosaic=0.01,
		   df0=50,
		   prob.nonMendelian=1.5e-6,
		   returnEmission=FALSE){ ## ignored
	stopifnot(states == 0:4)
	stopifnot(length(unique(ranges$chrom)) == 1)
	##family.id <- unique(ss(ranges$id))
	family.id <- unique(ranges$id)
	##fmonames <- paste(ss(family.id), c("03", "02", "01"), sep="_")
	pd2 <- phenoData2(trioSet)
	i <- match(family.id, sampleNames(trioSet))
	j <- match("CIDR_Name", colnames(pd2))
	stopifnot(!missing(i) && !missing(j))
	fmonames <- pd2[i, j, ]
	object <- computeLoglik(id=fmonames,
				trioSet=trioSet,
				ranges=ranges,
				mu.logr=mu.logr,
				states=states,
				baf.sds=baf.sds,
				prOutlier.logR=prOutlier[1],
				prOutlier.baf=prOutlier[2],
				prMosaic=prMosaic,
				df0=df0,
				returnEmission=returnEmission)
	if(returnEmission) return(object)
	trio.states <- trioStates(states)
	##tmp <- matrix(NA, nrow(trio.states), 2)
	tmp <- rep(NA, nrow(trio.states))
	##colnames(tmp) <- c("DN=0", "DN=1")
	## take into account 'the prior'
	##  -- the probability of transitioning to and from an altered state for each range
	##  -- the initial state probability if range is denovo
	##  -- the initial state probability of range is not denovo
	state.prev <- NULL
	denovo.prev <- NULL
	table1 <- readTable1(a=a)
	table3 <- readTable3(a=a)
	weightR <- 1/2
	state.names <- trioStateNames()
	norm.index <- which(state.names=="333")
	ranges <- ranges[order(start(ranges)), ]
	for(i in seq(length=nrow(ranges))){
		##if(i==7) break()
		obj <- object[which(range.index(object) == i), ]
		if(nrow(obj) < 2){
			##state.prev <- trio.states[norm.index, ]
			next()
		}
		LLR <- loglik(obj)["logR", , ,  ]
		LLB <- loglik(obj)["baf", , , ]
		LL <- weightR * LLR + (1-weightR)*LLB
		LLT <- matrix(NA, 3, 5)
		for(j in 1:3) LLT[j, ] <- apply(LL[, j, ], 2, sum, na.rm=TRUE)
		rownames(LLT) <- c("F", "M", "O")
		colnames(LLT) <- paste("CN_", states, sep="")
		if(FALSE){
			ii <- which(range.index(object)==i)
			apply(LLB[, 1, ], 2, sum, na.rm=TRUE)
			apply(LLB[, 1, ], 2, sum, na.rm=TRUE)
			apply(LLB[, 3, ], 2, sum, na.rm=TRUE)
			apply(LLR[, 3, ], 2, sum, na.rm=TRUE)
			r <- logR(obj)[, 3]
			tmp <- cbind(LLR[, 3, c(2,3)], r)
		}
		for(j in 1:nrow(trio.states)){
			tmp[j] <- joint1(LLT=LLT,
					 trio.states=trio.states,
					 tau=tau,
					 log.pi=log.pi,
					 normal.index=normal.index,
					 segment.index=i,
					 state.index=j,
					 table1=table1,
					 table3=table3,
					 state.prev=state.prev,
					 prob.nonMendelian=prob.nonMendelian)
		}
		## RS 4/29/2011
		argmax <- which.max(tmp)
		lik.norm <- tmp[norm.index]
		ranges$lik.norm[i] <- lik.norm
		ranges$argmax[i] <- argmax
		ranges$lik.state[i] <- tmp[argmax]
		stopifnot(!is.na(ranges$lik.state[i]))
		state.prev <- trio.states[argmax, ]
	}
	ranges$state <- trioStateNames()[ranges$argmax]
	ranges
}

calculateMarkerSd <- function(object, marker.index, sample.index){
	## exclude wga samples
	open(logR(object))
	##autosome.index <- which(chromosome(object) %in% 1:22)
	chr <- chromosome(object)[autosome.index]
	marker.index <- splitIndicesByLength(marker.index, ocProbesets())
	##j <- which(!sssampleNames(object) %in% wga.family)
	j <- sample.index
	message("excluding ", ncol(bsSet)-length(j), " WGA/CIDR samples")
	sds <- ocLapply(marker.index, function(i, object, j){
		lr <- as.matrix(logR(object)[i, j])
		sds <- rowMAD(lr, na.rm=TRUE)
		return(sds)
	}, object=object, j=j)
	close(logR(object))
	sds <- unlist(sds)
	##i <- unlist(marker.index)
	##fData(object)$MAD[i] <- sds
	return(sds)
	##object
}

pennTable <- function(a=0.0009, M=array(as.double(0), dim=rep(5,6))){
	res <- .C("calculateCHIT", a=as.double(a), M=M)
}


gridsetup <- function(figname, lattice.object, rd, ...){
	if(!missing(figname))
		trellis.device(device="pdf", file=figname, onefile=FALSE, ...)
	stopifnot(!missing(rd))
	chr.name <- paste("chr", rd$chrom[[1]], sep="")
	for(i in seq_along(lattice.object)){
		grid.newpage()
		lvp <- viewport(x=0, width=unit(0.52, "npc"), just="left", name="lvp")
		pushViewport(lvp)
		pushViewport(dataViewport(xscale=lattice.object[[i]][[1]]$x.limits,
					  yscale=c(0,1), clip="on"))
		print(lattice.object[[i]][[1]], newpage=FALSE, prefix="plot1", more=TRUE)
		upViewport(0)
		grid.text("Log R Ratio", x=unit(0.25, "npc"), y=unit(0.96, "npc"), gp=gpar("cex"=0.8))
		grid.text("B allele frequency", x=unit(0.75, "npc"), y=unit(0.96, "npc"), gp=gpar("cex"=0.8))
		grid.text(paste(chr.name, ", Family", ss(rd$id[i])), x=unit(0.5, "npc"), y=unit(0.98, "npc"), gp=gpar("cex"=0.9))
		upViewport(0)
		print(lattice.object[[i]][[2]], position=c(0.48, 0, 1, 1), more=TRUE, prefix="baf")
	}
	if(!missing(figname)) dev.off()
	TRUE
}

## rewrite this in C
joint1c <- function(loglik,##object,
		   trio.states,
		   tau,
		   log.pi,
		   normal.index,
		   segment.index,
		   state.index,
		   table1,
		   table3,
		   is.denovo,
		   Prob.DN=1.5e-6,
		   denovo.prev,
		   state.prev) {
	state <- trio.states[state.index, ]
	##fmo <- list()
	fmo <- matrix(NA, nrow(loglik), 3)
	##tmp <- as.matrix(do.call("cbind", fmo))
	##for(i in 1:3) fmo[, i] <- loglik(object)["logR", , i, state[i]] + loglik(object)["baf", , i, state[i]]
	for(i in 1:3) fmo[, i] <- loglik["logR", , i, state[i]] + loglik["baf", , i, state[i]]
	fmo <- fmo[rowSums(is.na(fmo)) == 0, ]
	##fmo <- fmo[[1]]+fmo[[2]]+fmo[[3]]
	##f.m.o <- sum(fmo, na.rm=TRUE)
	if(segment.index == 1){
		## assume Pr(z_1,f | lambda) = Pr(z_2,m | lambda) = pi
		## For offspring, we have Pr(z_1,o | z_1,f, z_1,m, DN=0, 1)
		##    or 1/5 if DN=1
		##
		## if DN is 0 (not devovo), then many of the hidden
		##  states should have essentially an epsilon
		##  probability of occurring.
		pi.offspring <- c(lookUpTable1(table1, state),  1/5)
		pi.offspring <- pi.offspring[[is.denovo+1]]
		## could log.pi just have length(0:4)??
		log.pi <- c(log.pi[state[1]], log.pi[state[2]], log(pi.offspring))
		fmo <- apply(fmo, 2, sum, na.rm=TRUE)
		fmo <- fmo + log.pi
##		for(j in 1:2) fmo[j] <- fmo[j]+log.pi[state.index]
##		f <- sum(fmo[[1]])+log.pi[state.index]
##		m <- sum(fmo[[2]])+log.pi[state.index]
##		o <- sum(fmo[[3]])+log(pi.offspring)
		## p.E138 top left: Pr(DN_1 = 1 | model)
		##  -- the probability the first marker is in a denovo region
		log.Prob.DN <- ifelse(is.denovo, log(Prob.DN), log(1-Prob.DN))
		fmo[3] <- fmo[3]+log.Prob.DN
	} else{
		##** Note that when state is the 'normal.state'**
		##    (state.index == normal.index)
		##    tau is the probability of staying in the same state
		##
		##for k = normal.index, it would do the right thing
		## prob. leaving normal state to state k
		fmo <- apply(fmo, 2, sum, na.rm=TRUE)
		for(j in 1:2) fmo[j] <- fmo[j]+log(tau[state.prev[j], state[j]])

		##f <- log(tau[state.prev[1], state[1]]) + sum(fmo[[1]])
		##m <- log(tau[state.prev[2], state[1]]) + sum(fmo[[2]])
		if(denovo.prev & is.denovo){
			## Equation 9: Wang et al.
			## if previous marker was in a range that was called denovo
			## Pr(z_j,o, z_j-1,o | DN_j=1, DN_j-1=1)
			##  = Pr(z_j,o | z_j-1,o, DN_j=1, DN_j-1=1) * Pr(z_j-1,o|DN_j-1=1)
			##    2nd term is 1/5
			##    1st term is
			## Pr(CN state of offspring in previous segment | previous segment was denovo) = 1/5
			fmo[3] <- log(1/5) + log(tau[state.prev[3], state[3]]) + fmo[3]
		}
		if(denovo.prev & !is.denovo){
			## Equation 10: Wang et al.
			## previous marker was in a range that was not called denovo
			tabled.value <- lookUpTable1(table1, state)
			fmo[3] <- log(1/5) + log(tabled.value) + fmo[3]
		}
		if(!denovo.prev & is.denovo){
			## Equation 10
			tabled.value <- lookUpTable1(table1, state.prev)
			fmo[3] <- log(1/5) + log(tabled.value) + fmo[3]
		}
		if(!denovo.prev & !is.denovo){
			## 25 x 25 x 25 table available in source code of software
			tabled.value <- lookUpTable3(table3, state.prev, state.curr=state)
			fmo[3] <- log(tabled.value) + fmo[3]
		}
	}
	## add the probability of transitioning back to the normal state
	for(j in 1:3) fmo[j] <- fmo[j] + log(tau[state[j], 3])
	##f <- f+log(tau[state[1], 3])
	##m <- m+log(tau[state[2], 3])
	##o <- o+log(tau[state[3], 3])
	res <- sum(fmo)
	stopifnot(!is.na(res))
	##
	## we all need a transition proability for not denovo -> denovo -> not denovo
	## (Equation 11)
	##f.m.o <- sum(fmo, na.rm=TRUE)
	## prob. back to  normal state
	##f.m.o <- f.m.o+tau[state.index, normal.index]
	return(res)
}

completeTrios <- function(bsSet){
	is.father <- which(bsSet$father != "0")
	is.mother <- which(bsSet$mother != "0")
	dup.index <- grep("240", sssampleNames(bsSet)[is.father])
	is.father <- is.father[-dup.index]
	is.mother <- is.mother[-dup.index]
	family <- sssampleNames(bsSet)[is.father]
	trio.matrix <- cbind(bsSet$father[is.father],
			     bsSet$mother[is.mother],
			     ssampleNames(bsSet)[is.father])
	colnames(trio.matrix) <- c("F", "M", "O")
	rownames(trio.matrix) <- ss(trio.matrix[,1])
	##
	## we do not have data on all of the samples
	##
	nas <- match(as.character(trio.matrix), ssampleNames(bsSet))
	namatrix <- matrix(nas, nc=3)
	i <- which(rowSums(is.na(namatrix)) > 0)
	trio.matrix <- trio.matrix[-i, ]
	return(trio.matrix)
}

trioNames <- function(bsSet, pass.qc=TRUE){
	if(pass.qc){
		unique(sssampleNames(bsSet)[bsSet$complete.trio & bsSet$pass.qc & !bsSet$is.duplicate])
	} else{
		unique(sssampleNames(bsSet)[bsSet$complete.trio & !bsSet$is.duplicate])
	}
}

qcFlag <- function(bsSet){
	oligoClasses:::open(bsSet$MAD)
	qcFlag <- bsSet$is.wga | bsSet$MAD[] >= 0.3
	oligoClasses:::close(bsSet$MAD)
	sns.flag <- unique(sssampleNames(bsSet)[qcFlag])
}

calculateDenovoFrequency <- function(ranges.md, penn.offspring, bychrom=FALSE){
	ranges.md$state <- gsub("5", "4", ranges.md$state)
	not.normal <- isDenovo(ranges.md$state)
	if(!any(not.normal)) {
		df <- NULL
	} else {
		df <- as.data.frame(table(ranges.md$state[not.normal]))
		df$method <- "mindist"
		colnames(df) <- c("state", "freq", "method")
	}
	stopifnot(all(penn.offspring$pass.qc))
	not.normal <- coverage(penn.offspring) >= 10 & isDenovo(penn.offspring$state)
	df2 <- as.data.frame(table(penn.offspring$state[not.normal]))
	df2$method <- "PennCNV"
	if(!bychrom){
		colnames(df2) <- c("state", "freq", "method")
	}
	if(!is.null(df)){
		df <- rbind(df, df2)
	} else df <- df2
	##df$is.denovo <- df$call %in% c(deletionStates(), duplicationStates())
	##df$is.denovo <- isDenovo(df$state)
	##df <- df[isDenovo(df$call), ]
	df$col <- rep("white", nrow(df))
	##df$col[df$is.denovo] <- "blue"
	i1 <- df$method=="mindist" ##& df$is.denovo
	if(any(i1)){
		f1 <- df$freq[i1]
		names(f1) <- as.character(df$state[i1])
		f1 <- f1[order(names(f1))]
	} else f1 <- NULL
	i2 <- df$method=="PennCNV"## & df$is.denovo
	f2 <- df$freq[i2]
	names(f2) <- as.character(df$state[i2])
	f2 <- f2[order(names(f2))]
	if(!is.null(f1)){
		f1 <- f1[names(f1) %in% names(f2)]
		f2 <- f2[names(f2) %in% names(f1)]
		stopifnot(all.equal(names(f1), names(f2)))
		f1 <- f1[match(names(f1), names(f2))]
		df <- data.frame(mindist=f1, penn=f2)
		rownames(df) <- names(f1)
	} else {
		rownames(df) <- names(f2)
		df2 <- df
		df2$freq <- 0
		df2$method="mindist"
		df <- rbind(df2, df)
	}
	return(df)
}

matchingDiscordantRanges <- function(penn.offspring, ranges.md, top=1:100, state=332, min.coverage=10){
	##colnames(penn.offspring)
	penn.offspring$state <- harmonizeStates(penn.offspring, filter.multistate=TRUE)
	penn332 <- penn.offspring[penn.offspring$state==state & penn.offspring$pass.qc & penn.offspring$nmarkers >= min.coverage, ]
	penn332 <- penn332[order(penn332$nmarkers, decreasing=TRUE), ]
	penn332 <- penn332[top, ]
	mind <- ranges.md[ranges.md$argmax != 61, ]
	mdMatch <- vector("list", length(top))
	##mind33 <- ranges.md[ranges.md$state == 332 | ranges.md$state ==
	## exclude the ranges that overlap with mindist
	chrom <- unique(penn332$chrom)
	m <- 1
	for(i in 1:nrow(penn332)){
		query <- IRanges(start(penn332)[i], end(penn332)[i])
		md <- mind[mind$id == penn332$id[i] & mind$chrom == penn332$chrom[i], ]
		if(nrow(md) == 0) next()
		md$pennIndex <- i
		subj <- IRanges(start(md), end(md))
		mm <- matchMatrix(findOverlaps(query, subj))
		##stopifnot(nrow(mm) <= 1)
		if(nrow(mm) == 0) next()
		ii <- unique(as.integer(mm[, "subject"]))
		mdMatch[[i]] <- md[ii, ]
	}
	return(list(penn=penn332, md.matching=mdMatch))
}

argmax2state <- function(argmax){
	trio.states <- trioStates()
	trio.states <- paste(trio.states[,1], trio.states[,2], trio.states[,3], sep="")
	states <- trio.states[argmax]
	return(states)
}

discordance <- function(rd1, rd2, I.STATE, ...){
	##rd1 is query
	##rd2 is subject
	stopifnot(!missing(I.STATE))
	stopifnot(all(rd1$state==I.STATE))
	nOverlap <- matchingCall <- isConcordant <- rep(NA, nrow(rd1))
	chrom.id <- paste(rd1$chrom, rd1$id, sep="_")
	uid <- unique(chrom.id)
	for(i in seq_along(uid)){
		if(i %% 100==0) message(i, "\n")
		ix <- which(chrom.id==uid[i])
		qrd <- rd1[ix, ]
		chrom <- unique(qrd$chrom)
		id <- unique(qrd$id)
		stopifnot(length(id)==1 & length(chrom)==1)
		srd <- rd2[rd2$chrom==chrom & rd2$id == id, ]
		qr <- IRanges(start(qrd), end(qrd))
		sr <- IRanges(start(srd), end(srd))
		if(nrow(srd) == 0){
			browser()
			isConcordant[ix] <- matchingCall[ix] <- NA
			nOverlap[ix] <- 0
			next()
		}
		nn <- countOverlaps(qr, sr)
		nOverlap[ix] <- nn
		if(all(nn==0)){
			next()
		}
		ix <- ix[nn > 0]
		## now, find the matching ranges for ranges with at least 1 overlap
		mm <- matchMatrix(findOverlaps(qr, sr))
		##sapply(split(seq_along(mm[,1]), mm[,1]), length)
		##nOverlap[ix] <- nrow(mm)
		stopifnot(nrow(mm) > 0)
		state <- srd$state[mm[, "subject"]]
		state.spl <- split(state, mm[, "query"])
		state.cat <- sapply(state.spl, function(x) paste(x,collapse="-"))
		state.cat <- unlist(state.cat)
		matchingCall[ix] <- state.cat
		##isConcordant[ix] <- any(argmax == ARGMAX)
		isConc <- sapply(state.spl, function(x, I.STATE) any(x == I.STATE), I.STATE=I.STATE)
		isConcordant[ix] <- isConc
		if(i == length(uid)) cat("\n")
	}
	rd1$isConcordant <- isConcordant
	rd1$nOverlap <- nOverlap
	rd1$matchingCall <- matchingCall
	return(rd1)
}

initializeTrioContainer <- function(path, samplesheet, pedigree, trio.phenodata, chromosomes=1:22, cdfName, ..., verbose){
	stopifnot(require(ff))
	stopifnot(all(chromosomes %in% 1:22))
	##stopifnot(all(file.exists(dirname(filenames))))
	stopifnot(all(file.exists(file.path(path, paste(rownames(samplesheet), ".txt", sep="")))))
	stopifnot(!missing(pedigree))
	if(is.null(rownames(pedigree))){
		rns <- apply(pedigree, 1, paste, collapse=",")
		rownames(pedigree) <- rns
	}
	stopifnot(!any(duplicated(rownames(pedigree))))
	##samplesheet <- samplesheet[-match(c("sampleNames", "filenames"), names(samplesheet))]
	fD <- constructFeatureData(list.files(path, full.names=TRUE)[1],
					   cdfName=cdfName)
	ss <- array(NA, dim=c(nrow(pedigree), ncol(samplesheet), 3),
		    dimnames=list(rownames(pedigree),
		    colnames(samplesheet),
		    c("F", "M", "O")))
	father.index <- match(pedigree[, "F"], s(samplesheet$Sample.Name))
	mother.index <- match(pedigree[, "M"], s(samplesheet$Sample.Name))
	offspring.index <- match(pedigree[, "O"], s(samplesheet$Sample.Name))
	ss[, , "F"] <- as.matrix(samplesheet[father.index, ])
	ss[, , "M"] <- as.matrix(samplesheet[mother.index, ])
	ss[, , "O"] <- as.matrix(samplesheet[offspring.index, ])
	marker.index.list <- split(seq(length=nrow(fD)), fD$chromosome)
	stopifnot(all(diff(order(fD$chromosome, fD$position))>0))
	trioSets <- vector("list", length(chromosomes))

	for(j in seq_along(chromosomes)){
		CHR <- chromosomes[j]
		if(verbose) message("\t Chromosome ", CHR)
		L <- length(marker.index.list[[CHR]])
		logR <- createFF(paste("logR_chr", CHR, "_", sep=""),
				 dim=c(L, nrow(pedigree), 3),
				 vmode="double")
		baf <- createFF(paste("baf_chr", CHR, "_", sep=""),
				dim=c(L, nrow(pedigree), 3),
				vmode="double")
		fd <- fD[marker.index.list[[j]], ]
		dimnames(logR) <- dimnames(baf) <- list(featureNames(fd),
							as.character(pedigree[, "O"]),
							c("F", "M", "O"))
		pD <- annotatedDataFrameFrom(logR[, , 1], byrow=FALSE)
		trioSets[[j]] <- new("TrioSet", logRRatio=logR,
				     BAF=baf,
				     phenoData=pD,
				     featureData=fD[marker.index.list[[CHR]], ],
				     mindist=NULL,
				     annotation=cdfName)
		##sampleNames(phenoData(trioSets[[CHR]])) <- pedigree[, "O"]
		## add data to phenoData2 slot
		## (note: parents with multiple sibs are repeated)
		trioSets[[CHR]]@phenoData2 <- ss
	}
	return(trioSets)
}

minimumDistance <- function(path, samplesheet, pedigree,
			    container.filename,
			    chromosomes=1:22,
			    cdfName, file.ext=".txt",
			    readFiles=TRUE,
			    calculate.md=TRUE,
			    calculate.mad=TRUE,
			    exclusionRule, ## for calculateing row-wise mads
			    ..., ##additional arguments for segment
			    verbose=TRUE){#samplesheet, ...){
	stopifnot(nrow(pedigree) > 1) ## need to fix initialization of trioSet object otherwise
	## the rownames of the samplesheet correspond to the name of the parsed beadstudio data
	stopifnot(all(file.exists(file.path(path, paste(rownames(samplesheet), file.ext, sep="")))))
	if(!file.exists(container.filename)){
	##---------------------------------------------------------------------------
	##
	## initialize container
	##
	##----------------------------------------------------------------------------
		stopifnot(file.exists(dirname(container.filename)))
		stopifnot("Sample.Name" %in% colnames(samplesheet))
		stopifnot(colnames(pedigree) == c("F", "M", "O"))
		if(verbose) message("Instantiating a container for the assay data.")
		container <- initializeTrioContainer(path,
						     samplesheet,
						     pedigree,
						     trio.phenodata,
						     chromosomes, cdfName,
						     verbose=verbose)
		if(verbose) message("Saving as ", container.filename)
		container <- as(container, "TrioSetList")
		save(container, file=container.filename)
	} else {
		load(container.filename)
		container <- get("container")
	}
	##---------------------------------------------------------------------------
	##
	## reading processed files
	##
	##----------------------------------------------------------------------------
	if(readFiles){
		if(verbose) message("Reading ", nrow(pedigree), " files")
		##fatherIds <- as.character(paste(fullId(container[[1]])[, "F"], file.ext, sep=""))
		mads <- matrix(NA, nrow(pedigree), 3)
		dimnames(mads) <- list(sampleNames(container), c("F", "M", "O"))
		mads[, "F"] <- readParsedFiles(path, "F", container, chromosomes, file.ext, verbose)
		##motherIds <- as.character(paste(fullId(container[[1]])[, "M"], file.ext, sep=""))
		mads[, "M"] <- readParsedFiles(path, "M", container, chromosomes, file.ext, verbose)
		##offspringIds <- as.character(paste(fullId(container[[1]])[, "O"], file.ext, sep=""))
		mads[, "O"] <- readParsedFiles(path, "O", container, chromosomes, file.ext, verbose)
		mad(container[[1]]) <- mads
		for(CHR in 2:22){
			container[[CHR]]@mad <- mad(container[[1]])
		}
		save(container, file=container.filename)
	} else {
		if(verbose) message("readFiles is FALSE.")
	}
	##---------------------------------------------------------------------------
	##
	## calculate minimum distance
	##
	##----------------------------------------------------------------------------
	if(calculate.md){
		if(verbose) message("Computing the minimum distance")
		container <- lapply(container, calculateMindist)
		container <- as(container, "TrioSetList")
	}
	##---------------------------------------------------------------------------
	##
	## Calculate the MAD of the log R ratios for every sample
	##    -> ntrios x 3 matrix
	##    -> put in slot 'mad'
	##---------------------------------------------------------------------------
	if(calculate.mad){
		container <- calculateMads(container, exclusionRule, chromosomes, verbose)
		if(verbose) message("\tSaving updated container to ", container.filename)
		save(container, file=container.filename)
	}
	return(container)
}

minimumDistanceCalls <- function(id, container,
				 chromosomes=1:22,
				 ranges,
				 cbs.filename,
				 segment.md=missing(cbs.filename)&missing(ranges),
				 offspring.segs,
				 mindistance.threshold=0.075,
				 calculate.lr=TRUE,
				 prOutlier=c(0.01, 1e-15),
				 prMosaic=0.01,
				 prob.nonMendelian=1.5e-6,
				 mu.logr=c(-2, -0.5, 0, 0.3, 0.75),
				 ##baf.sds=c(0.02, 0.3, 0.02),
				 baf.sds=c(0.005, 0.2, 0.005),
				 pruneMinimumDistance=FALSE,
				 prune.by.call=TRUE,
				 min.coverage=10,
				 returnEmission=FALSE,
				 ..., verbose=TRUE, DNAcopy.verbose=0){
	##---------------------------------------------------------------------------
	##
	## Segment the minimum distance
	##
	##---------------------------------------------------------------------------
	if(missing(id)) {
		id <- sampleNames(container)
	} else stopifnot(all(id %in% sampleNames(container)))
	stopifnot(all(chromosomes %in% 1:22))
	if(missing(ranges)){
		if(!missing(cbs.filename)){
			if(file.exists(cbs.filename)) {
				loadRanges <- TRUE
				segment.md <- FALSE
			} else {
				loadRanges <- FALSE
				segment.md <- TRUE
			}
		} else {
			stop("Must provide name of file to store results from circular binary segmentation of the minimum distance")
		}
	} else {
		segment.md <- FALSE
		loadRanges <- FALSE
	}
##	if(missing(segment.md)){
##		if(file.exists(cbs.filename)) {
##			message("segment.md is missing but ", basename(cbs.filename), " already exists. Loading saved segmentation")
##			segment.md <- FALSE
##		} else{
##			segment.md <- TRUE
##		}
##	}
	if(segment.md){
		stopifnot(!missing(cbs.filename))
		stopifnot(file.exists(dirname(cbs.filename)))
		df <- xsegment(container[chromosomes], id=id, ..., verbose=verbose,
			       DNAcopy.verbose=DNAcopy.verbose)
		df$ID <- gsub("^[X]", "", df$ID)
		mdRanges <- RangedDataCBS(ranges=IRanges(df$loc.start, df$loc.end),
					  chromosome=df$chrom,
					  sampleId=df$ID,
					  coverage=df$num.mark,
					  seg.mean=df$seg.mean,
					  startIndexInChromosome=df$start.index,
					  endIndexInChromosome=df$end.index)
		if(!missing(offspring.segs)){
			message("\tChecking the offspring segmentation to see whether breakpoints occur within the minimum distance interval...")
			mdRanges <- narrow(mdRanges, offspring.segs, thr=mindistance.threshold)
			message("\tFinished 'narrowing' the minimum distance ranges")
		}
		mads <- container[[1]]$mindist.mad
		ix <- match(sampleNames(mdRanges), id)
		mdRanges$mindist.mad <- mads[ix]
		message("Saving the segmentation results from CBS (prior to pruning) to ", cbs.filename)
		save(mdRanges, file=cbs.filename)
	} else {
		if(loadRanges){
			if(missing(cbs.filename)) stop("cbs.filename is missing, but segment.md=FALSE")
			if(verbose) message("Loading saved cbs segmentation results")
			load(cbs.filename)
			nm <- strsplit(basename(cbs.filename), ".rda")[[1]][1]
			mdRanges <- get(nm)
		} else {
			mdRanges <- ranges
			rm(ranges);gc()
		}
	}
	##---------------------------------------------------------------------------
	##
	## Prune the minimum distance ranges
	##
	##---------------------------------------------------------------------------
	## compute likelihood ratio to infer most likely state
	if(calculate.lr){
		if(pruneMinimumDistance){
			message("Pruning ranges")
			pruned.segs <- prune(container[chromosomes],
					     ranges=mdRanges,
					     id=id,
					     lambda=0.05,
					     min.change=0.1,
					     min.coverage=10,
					     scale.exp=0.02,
					     verbose=verbose)
			rd <- stack(RangedDataList(pruned.segs))
			rd <- rd[, -grep("sample", colnames(rd))]
			prunedRanges <- RangedDataCBS(ranges=ranges(rd), values=values(rd))
			rm(rd, pruned.segs); gc()
		} else {
			message("Minimum distance ranges not pruned")
			prunedRanges <- mdRanges
			prunedRanges <- prunedRanges[sampleNames(prunedRanges) %in% id & chromosome(prunedRanges) %in% chromosomes, ]
			rm(mdRanges); gc()
		}
##		if(!missing(offspring.ranges)){
##			ii <- which(sampleNames(offspring.ranges) %in% id & chromosome(offspring.ranges) %in% chromosomes)
##			offspring.ranges <- offspring.ranges[ii, ]
##			if(verbose) message("Narrowing ranges")
##			prunedRanges <- narrow(prunedRanges, offspring.ranges)
##		}
		tau <- transitionProbability(states=0:4, epsilon=0.5)
		log.pi <- log(initialStateProbs(states=0:4, epsilon=0.5))
		message("Computing bayes factors")
		prunedRanges$state <- NA
		object <- computeBayesFactor(object=container[chromosomes],
					     ranges=prunedRanges,
					     tau=tau,
					     log.pi=log.pi,
					     prOutlier=prOutlier,
					     prMosaic=prMosaic,
					     prob.nonMendelian=prob.nonMendelian,
					     mu.logr=mu.logr,
					     baf.sds=baf.sds,
					     returnEmission=returnEmission)
		if(returnEmission){
			##return LogLik set
			return(object)
		} else {
			prunedRanges <- object
			rm(object);gc()
		}
		##---------------------------------------------------------------------------
		## The following line needs to be commented. Here's why
		##  suppose the x's indicate a CNV
		##   9----------xxxxxx--------30
		##  And the copy numbers are
		##   222222222220000002222222222
		##  If we remove small segments, we get
		##   22222222222 <gap>2222222222
		##  And after pruneByFactor, we would have
		##   222222222222222222222222222
		## prunedRanges <- prunedRanges[prunedRanges$num.mark >= min.coverage, ]
		## do a second round of pruning for adjacent segments
		## that have the same state
		if(prune.by.call){
			message("Pruning ranges")
			rd <- tryCatch(pruneByFactor(prunedRanges, f=prunedRanges$argmax, verbose=verbose),
				       error=function(e) NULL)
			if(is.null(rd)){
				message("Not able to pruneByFactor. Return the ranges without pruning")
				return(prunedRanges)
			}
			rd <- stack(RangedDataList(rd))
			prunedRanges <- rd[, -ncol(rd)]
		}
	}
	prunedRanges$state <- trioStateNames()[prunedRanges$argmax]
	prunedRanges <- prunedRanges[order(prunedRanges$id, prunedRanges$chrom, start(prunedRanges)), ]
	return(prunedRanges)
}

calculateMads <- function(container, exclusionRule, chromosomes, verbose){
	##---------------------------------------------------------------------------
	##
	## Calculate the MAD of minimum distance
	##    -> put in phenoData slot
	##
	##---------------------------------------------------------------------------
	message("Computing mad of the minimum distance.")
	sapply(container, function(x) invisible(open(mindist(x))))
	nc <- ncol(container[[1]])
	mads.md <- rep(NA, nc)
	for(j in 1:nc){
		m <- lapply(container, function(x, j) mindist(x)[, j], j=j)
		mm <- unlist(m)
		mads.md[j] <- mad(mm, na.rm=TRUE)
	}
	sapply(container, function(x) invisible(close(mindist(x))))
	## inefficient to put mad in each TrioSet if there is a large number of samples
	## just put in first
	if(verbose) message("\tStoring MAD in first element of the TrioSetList container")
	container[[1]]$mindist.mad <- mads.md
	rm(mads.md, m, mm); gc()
	## this is not a ff object, so we might want to update the .rda file here
	##---------------------------------------------------------------------------
	##
	## Calculate the MAD of the logR ratios (across samples) for each marker
	##    -> put in featureData slot
	##
	##---------------------------------------------------------------------------
	if(nc > 1){
		if(verbose) message("\tCalculating mad of log R ratios across offspring samples for each marker.  Can be helpful to exclude samples of low quality using the exclusionRule argument.")
		J <- seq(length=nc)
		if(!missing(exclusionRule)) {
			exclude.index <- exclusionRule(container[[1]])
			if(length(exclude.index) > 0){
				if(verbose) message("\t\tExcluding trio index ", exclude.index, " from MAD calculation")
				J <- J[-exclude.index]
			}
		}
		for(i in seq_along(container)){
			trioSet <- container[[i]]
			invisible(open(logR(trioSet)))
			lR <- logR(trioSet)[, J, "O"]
			if(is.matrix(lR)){
				if(ncol(lR) > 1){
					invisible(close(logR(trioSet)))
					fData(container[[i]])$marker.mad <- rowMAD(lR, na.rm=TRUE)
				}
			}
		}
	}
	return(container)
}

thresholdY <- function(object){
	ylim <- object$y.limits
	panel.args2 <- object$panel.args
	f <- function(args, ylim){
		y <- args$y
		if(any(y < ylim[1], na.rm=TRUE)){
			index <- which(y < ylim[1])
			n <- sum(y < ylim[1], na.rm=TRUE)
			y[index] <- ylim[1] + 0.2 + runif(n, -0.2, 0.2)
		}
		if(any(y > ylim[2], na.rm=TRUE)){
			index <- which(y > ylim[2])
			n <- sum(y > ylim[2], na.rm=TRUE)
			y[index] <- ylim[2] - 0.2 + runif(n, -0.2, 0.2)
		}
		args$y <- y
		return(args)
	}
	panel.args <- lapply(panel.args2, f, ylim)
}

xypanel <- function(x, y, panelLabels,
		    xlimit,
		    ylimit,
		    ##segments=TRUE,
		    ##segments.md=TRUE,
		    range, fmonames,
		    cbs.segs,
		    md.segs,
		    what,
		    ##ylim,  ylim is not passed
		    ..., subscripts){
	##panel.grid(v=10,h=10, "grey")
	panel.grid(v=0, h=4, "grey", lty=2)
	what <- unique(as.character(what[subscripts]))
	stopifnot(length(what) == 1)
	##if(panelLabels[panel.number()] == "min dist") y <- -1*y
	panel.xyplot(x, y, ...)
	index <- which(x >= start(range)/1e6 & x <= end(range)/1e6)
	##panelLabels <- rev(panelLabels)
	CHR <- range$chrom
	##pL <- as.character(panelLabels[panel.number()])
	if(what %in% c("father", "mother", "offspring")){
		if(!missing(cbs.segs)){
			segments <- TRUE
			id <- switch(what,
				     father=fmonames[1],
				     mother=fmonames[2],
				     offspring=fmonames[3])
			cbs.sub <- cbs.segs[sampleNames(cbs.segs)==id & chromosome(cbs.segs)==CHR, ]
		} else segments <- FALSE
	} else segments <- FALSE
	if(!missing(md.segs) & what == "min dist"){
		cbs.sub <- md.segs[md.segs$id %in% range$id, ]
		cbs.sub <- cbs.sub[cbs.sub$chrom == range$chrom, ]
		cbs.sub$seg.mean <- -1*cbs.sub$seg.mean
		segments.md <- TRUE
	} else segments.md <- FALSE
	if(segments | segments.md){
		if(missing(ylimit)) ylimit <- range(y, na.rm=TRUE) ##else ylim <- ylimit
		if(nrow(cbs.sub) > 0){
			index <- which(cbs.sub$seg.mean < ylimit[1])
			if(length(index) > 0)
				cbs.sub$seg.mean[index] <- ylimit[1] + 0.2
			index <- which(cbs.sub$seg.mean > ylimit[2])
			if(length(index) > 0)
				cbs.sub$seg.mean[index] <- ylimit[2] - 0.2
			stopifnot(nrow(cbs.sub) > 0)
			panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2, col="black")#gp=gpar("lwd"=2))
		}
	}
	if(what == "genes"){
		require(locuszoom)
		data(rf, package="locuszoom")
		rf <- rf[!duplicated(rf$geneName), ]
		rf.chr <- rf[rf$txStart/1e6 <= xlimit[2] & rf$txEnd/1e6 >= xlimit[1] & rf$chrom==paste("chr", CHR, sep=""), ]
		flatBed <- flatten.bed(rf.chr)
		flatBed$start <- flatBed$start/1e3
		flatBed$stop <- flatBed$stop/1e3
		panel.flatbed(flat=flatBed,
			      showIso=FALSE,
			      rows=5,
			      cex=0.6)
	}
	if(what=="CNV"){
		require(locuszoom)
		data(cnv, package="locuszoom")
		cnv.chr <- cnv[cnv$txStart/1e6 <= xlimit[2] & cnv$txEnd/1e6 >= xlimit[1] & cnv$chrom==paste("chr", CHR, sep=""), ]
		##cnv.chr$txStart=cnv.chr$txStart/1000
		##cnv.chr$txEnd=cnv.chr$txEnd/1000
		##current.viewport$xscale <- xlimit
		flatBed <- flatten.bed(cnv.chr)
		flatBed$start <- flatBed$start/1e3
		flatBed$stop <- flatBed$stop/1e3
		panel.flatbed(flat=flatBed,
			      showIso=FALSE, rows=5,
			      cex=0.6,
			      col="red")
	}
}

gridlayout <- function(figname, lattice.object, rd, cex.pch=0.3, ...){
	if(!missing(figname))
		trellis.device(device="pdf", file=figname, onefile=FALSE,
			       width=8, height=5)
	stopifnot(!missing(rd))
	chr.name <- paste("chr", rd$chrom[[1]], sep="")
	grid.newpage()
	lvp <- viewport(x=0, width=unit(0.50, "npc"), just="left", name="lvp")
	pushViewport(lvp)
	pushViewport(dataViewport(xscale=lattice.object[[1]]$x.limits,
				  yscale=c(0,1), clip="on"))
	print(lattice.object[[1]], newpage=FALSE, prefix="plot1", more=TRUE)
	## highlight points in range
	L <- seq(length=length(lattice.object[[1]]$panel.args))
	viewportNames <- paste("plot1.panel.1.", L, ".off.vp", sep="")
	panelArgs <- lattice.object[[1]]$panel.args[[1]]
	x <- panelArgs$x
	index <- which(x >= start(rd)/1e6 & x <= end(rd)/1e6)
##n	for(i in seq_along(viewportNames)){
##		seekViewport(rev(viewportNames)[i])
##		y <- lattice.object[[1]]$panel.args[[i]]$y
##		grid.points(x[index], y[index], pch=21,
##			    gp=gpar(cex=cex.pch, fill="lightblue", alpha=0.5))
####			    gp=gpar(cex=0.6, fill="lightblue"))
##	}
	seekViewport("plot1.panel.1.1.off.vp")
	grid.move.to(unit(start(rd)/1e6, "native"),
		     unit(0, "npc"))
	seekViewport(paste("plot1.panel.1.", max(L),".off.vp", sep=""))
	grid.line.to(unit(start(rd)/1e6, "native"),
		     unit(1, "npc"), ##gp=gpar(col="purple", lty="dashed", lwd=1))
		     gp=gpar(...))
	seekViewport("plot1.panel.1.1.off.vp")
	grid.move.to(unit(end(rd)/1e6, "native"),
		     unit(0, "npc"))
	seekViewport(paste("plot1.panel.1.", max(L),".off.vp", sep=""))
	grid.line.to(unit(end(rd)/1e6, "native"),
		     unit(1, "npc"),
		     gp=gpar(...))
		     ##gp=gpar(col="red", lwd=1, lty="dashed", col="purple"))
	upViewport(0)
	grid.text(paste(chr.name, ", Family", ss(rd$id)), x=unit(0.5, "npc"), y=unit(0.98, "npc"),
		  gp=gpar(cex=0.9))
	grid.text("Position (Mb)", x=unit(0.5, "npc"), y=unit(0.05, "npc"),
		  gp=gpar(cex=0.8))
	upViewport(0)
	print(lattice.object[[2]], position=c(0.5, 0, 0.98, 1), more=TRUE, prefix="plot2")
	L <- seq_along(lattice.object[[2]]$panel.args)
	seekViewport("plot2.panel.1.1.off.vp")
	grid.move.to(unit(start(rd)/1e6, "native"),
		     unit(0, "npc"))
	viewportName <- paste("plot2.panel.1.", max(L), ".off.vp", sep="")
	seekViewport(viewportName)
	grid.line.to(unit(start(rd)/1e6, "native"),
		     unit(1, "npc"),
		     gp=gpar(...))
	##gp=gpar(col="purple", lty="dashed", lwd=2))
	seekViewport("plot2.panel.1.1.off.vp")
	grid.move.to(unit(end(rd)/1e6, "native"),
		     unit(0, "npc"))
	seekViewport(viewportName)
	grid.line.to(unit(end(rd)/1e6, "native"),
		     unit(1, "npc"), ##gp=gpar(col="red", lwd=2, lty="dashed", col="purple"))
		     gp=gpar(...))
	L <- seq(length=length(lattice.object[[2]]$panel.args))
	viewportNames <- paste("plot2.panel.1.", L, ".off.vp", sep="")
	panelArgs <- lattice.object[[2]]$panel.args
	x <- panelArgs[[1]]$x
	index <- which(x >= start(rd)/1e6 & x <= end(rd)/1e6)
##	for(i in seq_along(viewportNames)){
##		seekViewport(rev(viewportNames)[i])
##		y <- panelArgs[[i]]$y
##		grid.points(x[index], y[index], pch=21,
##			    gp=gpar(cex=cex.pch, fill="lightblue", alpha=0.5))
##	}
	if(!missing(figname)) dev.off()
	TRUE
}

gridlayout2 <- function(method1, xyList, otherCall, ranges, cex.call=0.6){
	stopifnot(!missing(method1))
	stopifnot(method1 %in% c("penn", "md"))
	f1 <- xyList[[1]]
	f2 <- xyList[[2]]
	nms <- names(f1)
	nms2 <- paste(sampleNames(ranges), start(ranges), sep="_")
	ix <- match(nms, nms2)
	stopifnot(length(ix) > 0)
	stopifnot(all(nms2[ix] == nms))
	ranges <- ranges[ix, ]
	if(is(otherCall, "character")){
		if(otherCall == "no overlap")
			plotOtherCall <- FALSE
	} else plotOtherCall <- TRUE
	for(i in seq_along(f1)){
		gridlayout(lattice.object=list(f1[[i]], f2[[i]]),
			   rd=ranges[i, ],
			   col="orange", cex=0.5,
			   cex.pch=0.1, fill="transparent", lty="solid", lwd=1)
		if(is(otherCall, "character")){
			if(otherCall == "no overlap"){
				upViewport(0)
				grid.text("PennCNV call: 333",
					  x=unit(0.6, "npc"),
					  y=unit(0.01, "npc"),
					  gp=gpar(col="grey30", cex=cex.call))
				pr <- NULL
			}
		} else {
			pr <- otherCall[otherCall$id == ranges$id[i], ]
			if(nrow(pr) == 0){
				stopifnot(otherCall$method[1]=="PennCNV")
				upViewport(0)
				grid.text("PennCNV call: 333",
					  x=unit(0.6, "npc"),
					  y=unit(0.01, "npc"),
					  gp=gpar(col="blue", cex=cex.call))
				pr <- NULL
			}
		}
		if(!is.null(pr)){
			method <- unique(otherCall$method)
			if(method=="PennCNV"){
				pc <- paste(pr$triostate, collapse="-")
			}
			pc <- paste(state(pr), collapse="-")
			upViewport(0)
			grid.text(paste(method, "call:", pc),
				  x=unit(0.6, "npc"),
				  y=unit(0.01, "npc"),
				  gp=gpar(col="blue", cex=cex.call))
		}
		if(method1=="penn"){
			grid.text(paste(method1, "call:", ranges$triostate[i]),
				  x=unit(0.2, "npc"),
				  y=unit(0.01, "npc"),
				  gp=gpar(col="orange", cex=cex.call))
		}
		if(method1=="md"){
			grid.text(paste("min dist call:", state(ranges)[i]),
				  x=unit(0.2, "npc"),
				  y=unit(0.01, "npc"),
				  gp=gpar(col="orange", cex=cex.call))
		}
		lr1 <- ranges$lik.state[i]
		lr2 <- ranges$lik.norm[i]
		lrr <- lr1-lr2
		l <- formatC(c(lr1, lr2, lrr), digits=1, format="f")
		grid.text(paste("log(Pr(332 | .))=", l[1], "\n",
				"log(Pr(333 | .))=", l[2], "\n",
				"log ratio:", l[3], "\n"),
			  x=unit(1, "npc"),
			  y=unit(0, "npc"),
			  just=c("right", "bottom"),
			  gp=gpar(col="grey10", cex=cex.call))
		if(!is.null(pr)){
			pr <- pr[state(pr) != "333", ]
			if(nrow(pr)==0) next()
			L <- length(f2[[i]]$panel.args)
			upViewport(0)
			seekViewport("plot2.panel.1.1.off.vp")
			grid.segments(x0=unit(start(pr)/1e6, "native"),
				      y0=unit(-0.1, "npc"),
				      x1=unit(start(pr)/1e6, "native"),
				      y1=unit(0.05, "npc"),
				      gp=gpar(col="blue", lwd=2))
			grid.segments(x0=unit(end(pr)/1e6, "native"),
				      y0=unit(-0.1, "npc"),
				      x1=unit(end(pr)/1e6, "native"),
				      y1=unit(0.05, "npc"),
				      gp=gpar(col="blue", lwd=2))
			seekViewport("plot1.panel.1.1.off.vp")
			grid.segments(x0=unit(start(pr)/1e6, "native"),
				      y0=unit(-0.1, "npc"),
				      x1=unit(start(pr)/1e6, "native"),
				      y1=unit(0.05, "npc"),
				      gp=gpar(col="blue", lwd=2))
			grid.segments(x0=unit(end(pr)/1e6, "native"),
				      y0=unit(-0.1, "npc"),
				      x1=unit(end(pr)/1e6, "native"),
				      y1=unit(0.05, "npc"),
				      gp=gpar(col="blue", lwd=2))
		}
	}
}

splitByDistance <- function(x, thr=90e3){
	##factorList <- list()
	##for(i in 1:22){
        ##		x <- position(trioSets[[i]])
	f <- c(0, cumsum(diff(x) > thr))
	tab.f <- table(f)
	while(any(tab.f < 1000)){
		j <- which(tab.f < 1000)[[1]]
		factor.val <- as.integer(names(tab.f)[j])
		if(factor.val < max(f)){
			f[f==factor.val] <- factor.val+1
		} else {
			f[f==factor.val] <- factor.val-1
		}
		tab.f <- table(f)
	}
	return(f)
}

myfilter <- function(x, filter, ...){
	res <- filter(x, filter,...)
	##now fill out the NAs
	M <- length(filter)
	N <- (M- 1)/2
	L <- length(x)
	for(i in 1:N){
		w <- filter[(N-i+2):M]
		y <- x[1:(M-N+i-1)]
		res[i] <- sum(w*y)/sum(w)
		w <- rev(w)
		ii <- (L-(i-1))
		y <- x[(ii-N):L]
		res[ii] <- sum(w*y)/sum(w)
	}
	return(res)
}

minimumDistancePlot <- function(trioSets, ranges, md.segs, cbs.segs, frame=2e6){
	md.segs <- md.segs[chromosome(md.segs) %in% chromosome(ranges) & sampleNames(md.segs) %in% sampleNames(ranges), ]
	cbs.segs <- cbs.segs[chromosome(md.segs) %in% chromosome(ranges) & sssampleNames(md.segs) %in% sssampleNames(ranges),]
	r1 <- ranges
	f1 <- f2 <- list()
	for(i in 1:nrow(r1)){
		panelLabels <- c("father", "mother", "offspring", "min dist")
		f1[[i]] <- xyplot(r ~ x | id, trioSets, ##range=penn.offspring[index, ],
				  panel=xypanel,
				  range=r1[i, ],
				  panelLabels=panelLabels,
				  layout=c(1,length(panelLabels)),
				  index.cond=list(rev(seq_along(panelLabels))),
				  frame=frame,
				  md.segs=md.segs,
				  ylim=c(-2.5, 1),
				  par.settings=list(plot.symbol=list(pch=21, cex=0.2, col="grey50"),
				  par.xlab.text=list(cex=0.8),
				  par.ylab.text=list(cex=0.8)),
				  xlab="",
				  ylab="log R ratio",
				  scales=list(x=list(tick.number=10, cex=0.5, tck=c(1,0), axs="i"),
				  alternating=rep(1,3),
				  y=list(cex=0.5)),
				  par.strip.text=list(lines=0.8, cex=0.7),
				  cbs.segs=cbs.segs)
		f1[[i]]$panel.args <- MinimumDistance:::thresholdY(f1[[i]])
		##print(f1a)
		panelLabels <- c("father", "mother", "offspring")
		f2[[i]] <- xyplot(b ~ x | id, trioSets, range=r1[i, ],
				  panelLabels=panelLabels, frame=frame,
				  scales=list(x=list(cex=0.5, tick.number=10,
					      tck=c(1,0), alternating=1),
				  y=list(cex=0.5, tck=c(0,1),
				  alternating=c(2,2,2))),
				  ##alternating=c(0,0,2,2,2))),
				  layout=c(1,length(panelLabels)), index.cond=list(length(panelLabels):1), pch=21,
				  par.settings=list(plot.symbol=list(pch=21, cex=0.2, col="grey50"),
				  par.xlab.text=list(cex=0.8),
				  par.ylab.text=list(cex=0.8)),
				  par.strip.text=list(lines=0.8, cex=0.7),
				  xlim=f1[[i]]$x.limit,
				  xlab="", ylab="B Allele frequency")
	}
	names(f1) <- names(f2) <- paste(sampleNames(ranges), start(ranges), sep="_")
	return(list(f1, f2))
}

narrow <- function(md.range, cbs.segs, thr, verbose=TRUE){
	i <- which(sampleNames(cbs.segs) %in% sampleNames(md.range))
	cbs.segs <- cbs.segs[i, ]
	i <- which(chromosome(cbs.segs) %in% chromosome(md.range))
	if(length(i) > 0 & length(i) < nrow(cbs.segs)){
		cbs.segs <- cbs.segs[i, ]
	}
	chroms <- unique(chromosome(md.range))
	rdN <- list()
	if(verbose) {
		message("Narrowing the ranges with segment mean above ", thr, " if the offspring has a start/stop in the region")
		pb <- txtProgressBar(min=0, max=length(chroms), style=3)
	}
	for(i in seq_along(chroms)){
		if(verbose) setTxtProgressBar(pb, i)
		j <- chroms[i]
		md <- md.range[chromosome(md.range) == j, ]
		md <- md[order(sampleNames(md), start(md)), ]
		jj <- which(chromosome(cbs.segs)==j)
		cbs <- cbs.segs[jj, ]
		cbs <- cbs[order(sampleNames(cbs), start(cbs)), ]
		ir1 <- IRanges(start(md), end(md))
		ir2 <- IRanges(start(cbs), end(cbs))
		mm <- findOverlaps(ir1, ir2)
		qhits <- queryHits(mm)
		shits <- subjectHits(mm)
		index <- which(sampleNames(md)[qhits] == sampleNames(cbs)[shits])
		if(length(index) > 0){
			qhits <- qhits[index]
			shits <- shits[index]
		} else stop("no overlap")
		##---------------------------------------------------------------------------
		##
		## only narrow the range if the minimum distance is
		## bigger than some nominal value. Otherwise, we use
		## the minimum distance range as is.
		##
		##
		##---------------------------------------------------------------------------
		abs.thr <- abs(md$seg.mean) > thr
		## I1 is an indicator for whether to use the cbs start
		I1 <- start(cbs)[shits] >= start(md)[qhits] & start(cbs)[shits] <= end(md)[qhits] & abs.thr[qhits]
		I2 <- end(cbs)[shits] <= end(md)[qhits] & end(cbs)[shits] >= start(md)[qhits] & abs.thr[qhits]
		st <- start(cbs)[shits] * I1 + start(md)[qhits] * (1-I1)
		en <- end(cbs)[shits] * I2 + end(md)[qhits] * (1-I2)
		st.index <- (cbs$start.index[shits] * I1 + md$start.index[qhits]*(1-I1))
		en.index <- (cbs$end.index[shits] * I2 + md$end.index[qhits]*(1-I2))
		## For each md range, there should only be one I1 that is TRUE
		## If I1 and I2 are true, then a range is completely contained within the md segment
		ids <- md$id[qhits]
		##  |--------------|
		## ---|--------|-----
		## Becomes
		##  |-|--------|---|
		index <- which(I1 & I2)-1
		index <- index[index!=0]
		index <- index[ids[index] == ids[index+1]]
		if(length(index) > 1){
			en[index] <- st[index+1]-1
			en.index[index] <- st.index[index+1]-1
		}
		##split(I1, qhits)
		##split(I2, qhits)
		stopifnot(sapply(split(st, qhits), function(x) all(diff(x) >= 0)))
		nm <- apply(cbind(st.index, en.index), 1, function(x) length(x[1]:x[2]))
		## keep segment means the same as the minimum distance
		tmp <- RangedData(IRanges(st, en),
				    id=sampleNames(md)[qhits],
				    chrom=chromosome(md)[qhits],
				    num.mark=nm,
				    seg.mean=md$seg.mean[qhits],
				    start.index=st.index,
				    end.index=en.index,
				  mindist.mad=md$mindist.mad[qhits])
		##ranges.below.thr <- split(!abs.thr[qhits], qhits)
		##ns <- sapply(ranges.below.thr, sum)
		uid <- paste(tmp$id, start(tmp), tmp$chrom, sep="")
		duplicated(uid)
		tmp <- tmp[!duplicated(uid), ]
		## for each subject, the following must be true
		index <- which(tmp$id[-nrow(tmp)] == tmp$id[-1])
		stopifnot(all(end(tmp)[index] < start(tmp)[index+1]))
		rdN[[i]] <- tmp[order(tmp$id, start(tmp)), ]
	}
	if(verbose) close(pb)
	if(length(rdN) > 1){
		rdL <- stack(RangedDataList(rdN))
		rd <- rdL[, -grep("sample", colnames(rdL))]
		rdCbs <- RangedDataCBS(ranges=ranges(rd), values=values(rd))
	} else rdCbs <- RangedDataCBS(ranges=ranges(rdN[[1]]), values=values(rdN[[1]]))
	return(rdCbs)
}

plotRange <- function(range, trioSets, md.segs,cbs.segs, penn.offspring, frame=2e6){
	tmp <- MinimumDistance:::minimumDistancePlot(trioSets=trioSets,
						     ranges=range,
						     md.segs=md.segs,
						     cbs.segs=cbs.segs,
						     frame=frame)
	penn.call <- correspondingCall(range, penn.offspring, subject.method="PennCNV")
	MinimumDistance:::gridlayout2("md", tmp, penn.call, ranges=range)
}
