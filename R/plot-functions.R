sizeCategories <- function(nSnps){
	## group
	##(500, +], (200, 500], (100, 200], (50, 100], (25, 50], (10, 25], (5, 10], (3, 5], (0, 3]
	sizeCat <- rep("1 (0, 3]", length(nSnps))
	sizeCat[nSnps > 3 & nSnps <=5] <- "2 (3, 5]"
	sizeCat[nSnps > 5 & nSnps <=10] <- "3 (5, 10]"
	sizeCat[nSnps > 10 & nSnps <=25] <- "4 (10, 25]"
	sizeCat[nSnps > 25 & nSnps <=50] <- "5 (25, 50]"
	sizeCat[nSnps > 50 & nSnps <=100] <- "6 (50, 100]"
	sizeCat[nSnps > 100 & nSnps <=200] <- "7 (100, 200]"
	sizeCat[nSnps > 200 & nSnps <=500] <- "8 (200, 500]"
	sizeCat[nSnps > 500 & nSnps <= 1000] <- "9 (500, 1000]"
	sizeCat[nSnps > 1000] <- "10 (1000, +]"
	sizeCat
}

plotRegion <- function(x, outdir="/thumper/ctsa/beaty/scharpf/crlmmOut", windowsize=50,
		       allowHetParent=FALSE, cdfName="human610quadv1b",
		       ylim=c(0.5, 8)){
	data(samplesheet, package="Beaty")
	samplesheet$individualId <- getIndividualId2(samplesheet)
	platedir <- file.path(outdir, samplesheet[match(x[["sample"]], samplesheet$individualId), "Sample.Plate"])
	CHR <- x[["chr"]]
	fn <- file.path(platedir, paste("crlmmSetList_", CHR, ".rda", sep=""))
	load(file.path(fn))
	crlmmSetList <- get("crlmmSetList")
	if(length(crlmmSetList) == 3) crlmmSet <- as(crlmmSetList, "CrlmmSet")
	if(length(crlmmSetList) < 3){
		if(file.exists(file.path(platedir, paste("crlmmSet_", CHR, ".rda", sep="")))){
			message("loading crlmmSet")
			load(file.path(platedir, paste("crlmmSet_", CHR, ".rda", sep="")))
		} else {
			object <- as(crlmmSetList, "SnpCallSetPlus")
			if(!"isSnp" %in% fvarLabels(object)){
				fData(object)$isSnp <- isSnp(object)
			}
			message("computing copy number")
			cnOpts <- cnOptions(batch=rep("A", ncol(object)), cdfName=cdfName)
			object$batch <- cnOpts[["batch"]]
			##trace(oneBatch, browser)
			crlmmSet <- computeCopynumber(object, cnOptions=cnOpts)
			message("saving crlmmSet")
			save(crlmmSet, file=file.path(platedir, paste("crlmmSet_", CHR, ".rda", sep="")))
		}
	}
	rm(crlmmSetList); gc()
	fData(crlmmSet)$isSnp <- isSnp(crlmmSet)
	x$familyId <- substr(x[["sample"]], 1, 5)
	if(!"familyId" %in% varLabels(crlmmSet)){
		phenoData(crlmmSet) <- addSampleSheet(crlmmSet)
		pData(crlmmSet)$familyId <- getFamilyId(crlmmSet)
		pData(crlmmSet)$familyMember <- getIndividualId(crlmmSet)
	}
	trioSet <- crlmmSet[, crlmmSet$familyId==x[["familyId"]]]
	if(ncol(trioSet) > 3){
		trioSet <- trioSet[, trioSet$familyMember %in% c("01", "02", "03")]
	}
	if(ncol(trioSet) < 3) return(NULL)
	whoisit <- sapply(trioSet$familyMember, who)
	if(length(whoisit) != 3) stop("whoisit does not have length 3.  check familyMember variable in crlmmSet")
	trioSet <- trioSet[, match(c("father", "mother", "offspring"), whoisit)]
	sampleNames(trioSet) <- c("father", "mother", "offspring")
	whichIndices <- function(object, x){
		which(position(object) >= as.integer(x[["start"]]) & position(object) <= as.integer(x[["end"]]) & isSnp(object))
	}
	region <- index <- whichIndices(trioSet, x)
	if(length(region) < 1) return()
	gts <- calls(trioSet[region, ])
	gts <- cbind(gts, isBiparental.matrix(gts, allowHetParent=allowHetParent))
	colnames(gts)[4] <- "isBiparental"
	pHet <- apply(calls(trioSet[region, ]) == 2, 2, mean, na.rm=TRUE)
	pHom <- 1-pHet
	index <- (min(index)-windowsize):(max(index)+windowsize)
	index <- index[index >= 1 & index <= nrow(trioSet)]
	trioSet <- trioSet[index, ]
 	notBpiIndicator <- !isBiparental.SnpCallSetPlus(trioSet, allowHetParent=allowHetParent)
	y <- isBiparental.SnpCallSetPlus(trioSet, allowHetParent=allowHetParent)
	y <- y+1
	y[is.na(y)] <- 0
	yy <- y
	bpiOnly <- y == 2
	y <- jitter(y, amount=0.1)
	xx <- position(trioSet)
	plot(xx, y, pch=21, col=grey(.7), yaxt="n", xaxt="n", ylim=c(-.2, 2.2),
	     xlab="position (Mb)", ylab="", xaxt="n")
	points(xx[notBpiIndicator], y[notBpiIndicator], pch=21, bg="royalblue")
	axis(2, at=c(0,1,2), labels=c("NI", "notBPI", "BPI"))
	axis(1, at=pretty(xx), labels=pretty(xx/1e6))
	abline(v=c(x[["start"]], x[["end"]]), lty=2, col="royalblue")
	legend("topleft", bty="n", legend=x[["sample"]])
	mtext(paste("Chr", x[["chr"]]), side=3, outer=TRUE, line=0)
	CN <- copyNumber(trioSet)
	for(j in 1:ncol(CN)){
		plot(xx, CN[, j], pch=21, col=grey(.7), xaxt="n",
		     ylab="",
		     cex=0.7, log="y", ylim=ylim)
		points(xx[notBpiIndicator], CN[notBpiIndicator, j], pch=21, cex=0.8, bg="royalblue")
		abline(h=1:3, col=grey(0.6))
		legend("topleft", legend=c("father", "mother", "offspring")[j], bty="n")
		legend("topright", legend=paste("% Het:", round(pHet[j],2)), bty="n")
		legend("top", legend=paste("SNR:", round(trioSet$SNR[j],0)), bty="n")
		abline(v=c(x[["start"]], x[["end"]]), lty=2, col="royalblue")
	}
	gtConfs <- gtConfidence(trioSet)
	for(j in 1:3){
		if(j==1){
			##yy>0 means that the genotype was informative
			if(length(xx[yy>0]) > 0){
				plot(xx, gtConfs[, j], pch=".", xlab="position (Mb)", xaxt="n",
				     ylim=c(-0.02, 1.02), ylab="crlmm confidence")
				points(xx[yy>0], jitter(gtConfs[yy > 0, j], amount=0.02),
				       pch=21, cex=0.7, bg=c("black", "blue", "red")[j],
				       col=c("black", "blue", "red")[j])
			} else {
				plot(xx, rep(0, length(xx)), type="n", xlab="", ylab="")
			}
		} else {
			if(length(xx[yy>0]) > 0){
				points(xx[yy>0], jitter(gtConfs[yy>0, j], amount=0.02), pch=21, cex=0.7, bg=c("black", "blue", "red")[j],
				       col=c("black", "blue", "red")[j])
			}
		}
		if(j == 3) {
			abline(v=c(x[["start"]], x[["end"]]), lty=2, col="royalblue")
			legend("bottom", bty="n", col=c("black", "blue", "red"), pt.bg=c("black", "blue", "red"), pch=21, legend=c("F", "M", "O"))
		}
	}
	axis(1, at=pretty(xx), labels=pretty(xx/1e6), outer=TRUE)
	rm(crlmmSet); gc()
	return(gts)
}


myPlot <- function(rD, row, cnSet, surround, ylim=c(0.2,6), ...){
	##require(lattice)
	rangedData <- rD
	rD <- rD[row, ]
	start <- start(rD)
	end <- end(rD)
	index <- which(position(cnSet) >= start & position(cnSet) <= end)
	bpiIndex <- index
	f <- max(1, index[1] - surround)
	l <- min(nrow(cnSet), index[length(index)] + surround)
	index <- f:l
	##cnSet <- cnSet[index, ]
	notBpi <- !isBiparental.SnpSuperSet(cnSet[bpiIndex, ], allowHetParent=FALSE)
	notBpi <- notBpi==TRUE & !is.na(notBpi)
	y <- rev(c(0.75, 0.5, 0.25))
	x <- 1:sum(notBpi)
	plot(y=1:sum(notBpi), x=rep(y[1], length(x)),
	     type="n", xaxt="n", yaxt="n", xlab="", ylab="",
	     xlim=c(-0.2, 1.3))
	text(y=1:sum(notBpi), x=rep(y[1], length(x)),
	     labels=c("AA", "AB", "BB")[snpCall(cnSet)[bpiIndex[notBpi], 1]])
	for(j in 2:3){
		text(y=x, x=rep(y[j], length(x)), labels=c("AA", "AB", "BB")[snpCall(cnSet)[bpiIndex[notBpi], j]])
	}
	FMO <- snpCall(cnSet)[bpiIndex[notBpi], ]
	parentOfOrigin <- function(x){
		## F   M   O     Parent
		## AA  BB  AA    F
		## AB  AA  BB    F
		## AB  BB  AA    F
		## BB  AA  BB    F
		## AA  AB  BB    M
		## BB  AB  AA    M
		## AA  BB  BB    M
		## BB  AA  AA    M
		if(all(x == c(1, 3, 1)) |
		   all(x == c(2, 1, 3)) |
		   all(x == c(2, 3, 1)) |
		   all(x == c(3, 1, 3)))
			return("F")
		if(all(x == c(1, 2, 3)) |
		   all(x == c(3, 2, 1)) |
		   all(x == c(1, 3, 3)) |
		   all(x == c(3, 1, 1)))
			return("M")
		return("?")
	}
	transmittedFrom <- apply(FMO, 1, parentOfOrigin)
	text(y=x, x=rep(1, length(x)), labels=transmittedFrom, col="blue")
	axis(3, at=c(y,1), labels=c("F", "M", "O", "T"))
	cnSet <- cnSet[index, ]
	CN <- copyNumber(cnSet)
	y <- isBiparental.SnpSuperSet(cnSet, allowHetParent=FALSE)
	noty <- !y
	y <- y+1
	y[is.na(y)] <- 0
	yy <- y
	bpiOnly <- y == 2
	y <- jitter(y, amount=0.1)
	xx <- position(cnSet)
	plot(xx, y, pch=21, col=grey(.7), yaxt="n", xaxt="n", ylim=c(-.2, 2.2),
	     xlab="position (Mb)", ylab="", xaxt="n")
	points(xx[noty], y[noty], pch=21, bg="royalblue")
	axis(2, at=c(0,1,2), labels=c("NI", "notBPI", "BPI"))
	axis(1, at=pretty(xx), labels=pretty(xx/1e6))
	abline(v=c(start, end), lty=2, col="royalblue")
	legend("topleft", bty="n", legend=sampleNames(cnSet)[3])
	mtext(space(rD), side=3, outer=TRUE, line=0)
	CN <- copyNumber(cnSet)
	CN[CN < ylim[1]] <- ylim[1]
	CN[CN > ylim[2]] <- ylim[2]
	for(j in 1:ncol(CN)){
		plot(xx, CN[, j], pch=21, col=grey(.7), xaxt="n",
		     ylab="",
		     cex=0.7, ylim=ylim)
		points(xx[noty], CN[noty, j], pch=21, cex=0.8, bg="royalblue")
		abline(h=1:3, col=grey(0.6))
		legend("topleft", legend=c("father", "mother", "offspring")[j], bty="n")
		##legend("topright", legend=paste("% Het:", round(pHet[j],2)), bty="n")
		legend("top", legend=paste("SNR:", round(cnSet$SNR[j],1)), bty="n")
		abline(v=c(start, end), lty=2, col="royalblue")
	}
	gtConfs <- confs(cnSet)
	minConf <- apply(gtConfs, 1, min)
	minConf[minConf < 0.5] <- 0.5
	plot(xx, minConf, pch=21, col=grey(0.6), cex=0.8, xaxt="n", ylim=c(0.5, 1))
	axis(1, at=pretty(xx), labels=pretty(xx/1e6), outer=TRUE)
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
	legend("bottomright", legend=who(sampleName),  bty="n", cex=0.8)
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

plotSegs <- function(index, ranges1, ranges.md, ranges2, mset, FRAME, xlim, strict=FALSE, ylim=c(-1.5, 0.5), ...){
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
	ranges2 <- ranges2[family(ranges2) %in% family(this.range), ]
	ranges2$seg.mean[ranges2$seg.mean < ylim[1]] <- ylim[1]
	ranges2$seg.mean[ranges2$seg.mean > ylim[2]] <- ylim[2]
	## goal is to have the feature cover approx 5% of the plot
	if(missing(xlim)){
		if(missing(FRAME)){
			w <- width(this.range)
			FRAME <- w/0.05  * 1/2
		}
		i <- featuresInRange(mset, this.range, FRAME=FRAME)
	} else{
		i <- which(position(mset) >= xlim[1] & position(mset) <= xlim[2])
	}
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

getGenomicAxis <- function(deletion.ranges, CHR, FRAME=1e6, xlim){
	deletion.chr <- deletion.ranges[deletion.ranges$chrom == CHR, ]
	if(!any(deletion.chr$n.overlap > 0)){
		message("None of the segments on chromosome ", CHR, " have any overlap")
		return(integer(0))
	}
	nn <- deletion.chr$n.overlap
	index <- which(deletion.chr$n.overlap == max(nn))
	samples.chr <- unique(as.character(sapply(deletion.chr$others[index], function(x) strsplit(x, ", ")[[1]])))
	index.case <- deletion.chr$id[index]
	samples.chr <- unique(c(index.case, samples.chr))
	samples.chr <- samples.chr[samples.chr != "NA"]
	deletion.chr <- deletion.chr[deletion.chr$id %in% samples.chr, ]
	others <- unique(unlist(sapply(deletion.chr$others, function(x) strsplit(x, ", ")[[1]])))
	others <- others[others != "NA"]
	deletion.chr <- deletion.chr[deletion.chr$id %in% others, ]
	samples.chr <- unique(deletion.chr$id)
	## if samples has just one overlap, see if its within 200kb of the sample with most overlap
	position.chr <- c(min(start(deletion.chr)), max(end(deletion.chr)))
	if(missing(xlim)){
		xlim <- c(min(position.chr) -FRAME,
			  max(position.chr) +FRAME)
	}
	return(xlim)
}

plotCytobandWithRanges <- function(deletion.ranges, CHR, xlim, FRAME=1e6, ylab.colid, featureData, ...){
	deletion.22 <- deletion.ranges[deletion.ranges$chrom == CHR, ]
	nn <- deletion.22$n.overlap
	index <- which(deletion.22$n.overlap == max(nn))
	samples.22 <- unique(as.character(sapply(deletion.22$others[index], function(x) strsplit(x, ", ")[[1]])))
	index.case <- deletion.22$id[index]
	samples.22 <- unique(c(index.case, samples.22))
	samples.22 <- samples.22[samples.22 != "NA"]
	deletion.22 <- deletion.22[deletion.22$id %in% samples.22, ]
	others <- unique(unlist(sapply(deletion.22$others, function(x) strsplit(x, ", ")[[1]])))
	others <- others[others != "NA"]
	deletion.22 <- deletion.22[deletion.22$id %in% others, ]
	samples.22 <- unique(deletion.22$id)
	## if samples has just one overlap, see if its within 200kb of the sample with most overlap
	position.22 <- c(min(start(deletion.22)), max(end(deletion.22)))
	if(missing(xlim)){
		xlim <- c(min(position.22) -FRAME,
			  max(position.22) +FRAME)
	}
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
			if(nrow(this.range) == 1)
				text(xlim[2], mean(y), labels=this.range$num.mark[j], col="grey30", cex=0.6, adj=1)
		}
		if(nrow(this.range) > 1){
			text(xlim[2], mean(y), labels=median(this.range$num.mark), adj=1, cex=0.6, col="grey30")
		}
	}
	if(!missing(featureData)){
		fD <- featureData[featureData$chromosome == CHR, ]
		fD <- fD[fD$position >= xlim[1] & fD$position <= xlim[2]]
		rug(fD$position, side=3)
	}
	axis(2, at=seq_along(samples.22), labels=samples.22, cex.axis=0.6, adj=0)
	text(0, mean(position.22), paste(diff(position.22)/1e3, "kb"))
	return(list(xlim=xlim, v=position.22, samplesPlotted=samples.22))
}

myPanel <- function(x,y, subscripts, ...){
	panel.xyplot(x, y, ...)
	panel.grid(h = 10, v = 10, col = "grey", lty = 3)
}

rangesInXlim <- function(ranges.object, CHR, xlim){
	del <- ranges.object[ranges.object$chrom == CHR, ]
	del <- del[start(del) <= xlim[2] & end(del) >= xlim[1], ]
	return(del)
}


grid.rectangle <- function(ranges.object, CHR, xlim,
			   jugessur,
			   rf, cnv, cols, alpha=0.1,
			   deletion.colors=TRUE){
	chr.name <- paste("chr", CHR, sep="")
	del <- rangesInXlim(ranges.object, CHR, xlim)
	stopifnot(nrow(del) > 0)
	jug <- jugessur[start(jugessur) <= xlim[2] & end(jugessur) >= xlim[1] & jugessur$chrom==CHR, ]
	rf.chr <- rf[rf$txStart <= xlim[2] & rf$txEnd >= xlim[1] & rf$chrom==chr.name, ]
	cnv.chr <- cnv[cnv$txStart <= xlim[2] & cnv$txEnd >= xlim[1] & cnv$chrom==chr.name, ]
	flatBed <- flatten.bed(rf.chr)
	flatBed.cnv <- flatten.bed(cnv.chr)
	del <- del[order(del$id), ]
	ids <- del$id
	ids <- split(ids, ids)
	ids <- ids[match(unique(del$id), names(ids))]
	starts <- start(del)
	ends <- end(del)

	## for matching samples with another ranges.object
	## need to insert 'blank rows'
##	if(!missing(ranges.object2)){
##		uids <- unique(ranges.object2$id)
##		match(del$id, uids)
##		uidsNotInRanges <- uids[!uids %in% ranges.object$id]
##	}
	yy <- (1:(length(ids)+2))/(length(ids)+2)
	y <- yy[-c(1, length(yy))]
	y <- rep(y, sapply(ids, length))
	sp <- diff(yy[1:2])/2
	if(missing(cols)){
		message("color scheme for denovo homo vs hemi deletions")
		##cols <- brewer.pal(12, "Paired")[11:12]
		cols <- c("grey95", "grey55")
	}
	cols <- col2rgb(cols)/255
	cols <- rgb(cols[1, ],
		    cols[2, ],
		    cols[3, ], alpha=alpha)
	if(deletion.colors){
		homo.del <- substr(del$state, 3, 3) == 1  ##homozygous
		hemi.del <- substr(del$state, 3, 3) == 2  ##hemizygous
		fill <- rep(NA, nrow(del))
		fill[homo.del] <- cols[1]
		fill[hemi.del] <- cols[2]
	} else {
		denovo.amp <- substr(del$state, 3, 3) == 4
		fill[denovo.amp] <- cols[1]
		fill[!denovo.amp] <- cols[2]
	}
	## rectangle using graphical objects in grid
	layout <- grid.layout(3, 1, heights=unit(c(1, 0.2, 0.2), c("null", "null", "null")))
	##grid.newpage()
	pushViewport(viewport(layout=layout))
	pushViewport(viewport(layout.pos.row=1, name="row1"))
	pushViewport(plotViewport(c(2,4,2,2)))  ## rectangular region with space in lines to the bottom, left, top and right
	sts <- sapply(starts, function(x) x[[1]])
	eds <- sapply(ends, function(x) x[[1]])
	pushViewport(dataViewport(c(sts/1e6, unlist(ends)/1e6), rep(y, each=2), xscale=xlim/1e6,
				  yscale=c(0, 1), name="plotRegion"))
	grid.rect(x=starts/1e6, y=y, height=sp, width=width(del)/1e6,
		  default.units="native", gp=gpar(fill=fill, border=fill), name="cnv.regions", just=c("left", "center"))
	grid.rect(gp=gpar("col"="grey50"))
	grid.xaxis(gp=gpar("cex"=0.7), at=pretty(xlim/1e6, n=10))
	yax <- yaxisGrob(at=y, label=del$id, gp=gpar("cex"=0.7), name="ylabs")
	editGrob(yax, "labels", just="bottom")
	grid.draw(yax)
	panel.grid(h = -1, v = 10, col = "grey", lty = 3)
	upViewport()
	grid.text(paste("Chromosome", CHR), y=unit(1.05, "native"), just="top")
	if(nrow(jug) > 0){
		grid.rect(x=start(jug)/1e6, y=rep(yy[1], nrow(jug)), height=sp,
			  width=width(jug)/1e6,
			  default.units="native", gp=gpar(fill="lightblue"),
			  name="jug.regions", just=c("left", "top"), vp="plotRegion")
		grid.lines(y=yy[1]+sp, default.units="native", gp=gpar("lty"=2), vp="plotRegion")
		grid.text("Jugessur", gp=gpar("col"="blue"),
			  ##y=unit(yy[1]+sp, "native"),
			  y=unit(yy[1], "native"),
			  x=unit(xlim[1]/1e6, "native"), just=c("left", "bottom"), vp="plotRegion")
	}
	grid.text(label=as.character(del$num.mark), vp="plotRegion",
		  y=unit(y, "native"),
		  ##x=unit(ends/1e6+0.1, "native"),
		  x=unit((ends/1e6+starts/1e6)/2, "native"),
		  gp=gpar("cex"=0.7))
	upViewport(0)
	if(nrow(flatBed) > 0){
		mn <- min(flatBed$start)
		mx <- max(flatBed$stop)
		pushViewport(viewport(layout=layout))
		pushViewport(viewport(layout.pos.row=2, name="genes", clip="on"))
		pushViewport(plotViewport(c(0, 4, 0,2)))  ## needs to be exact same as above
		pushViewport(dataViewport(xscale=xlim/1000, yscale=c(0,1), clip="on"))
		panel.flatbed(flat=flatBed, showIso=FALSE, cex=0.6, rows=6)
	}
	if(nrow(flatBed.cnv) > 0){
		upViewport(0)
		mn <- min(flatBed.cnv$start)
		mx <- max(flatBed.cnv$stop)
		pushViewport(viewport(layout=layout))
		pushViewport(viewport(layout.pos.row=3, name="genes", clip="on"))
		pushViewport(plotViewport(c(0, 4, 0,2)))  ## needs to be exact same as above
		pushViewport(dataViewport(xscale=xlim/1000, yscale=c(0,1), clip="on"))
		panel.flatbed(flat=flatBed.cnv, showIso=FALSE, cex=0.6, rows=6, col="red", fill="red")
	}
	return()
}

grid.trioLogR <- function(mset,
			  labels=c("Father", "Mother", "Offspring", "max distance"),
			  cex=0.3, ...){
	grid.obj <- data.frame(logR=c(logR.F(mset),
			       logR.M(mset),
			       logR.O(mset),
			       -mindist(mset)),
			       x=rep(position(mset)/1e6, 4),
			       subject=factor(rep(1:4, each=nrow(mset)), labels=labels, ordered=TRUE))
	if("ylim" %in% names(list(...))){
		ylim <- list(...)[["ylim"]]
		nn <- sum(grid.obj$logR < ylim[1], na.rm=T)
		grid.obj$logR[grid.obj$logR < ylim[1]] <- runif(n=nn, min=ylim[1], max=ylim[1]+0.5)
	}
	nn <- sum(grid.obj$logR < -4, na.rm=TRUE)
	if(nn > 10) {
		ii <- which(grid.obj$logR < -5)
		grid.obj$logR[ii] <- jitter(rep(-5, length(ii)),  amount=0.5)
	} else {
		ii <- which(grid.obj$logR < -2)
		if(length(ii) > 0){
			grid.obj$logR[ii] <- jitter(rep(-2, length(ii)),  amount=0.5)
		}
	}
	nn <- sum(grid.obj$logR > 1, na.rm=TRUE)
	if(nn > 10) {
		ii <- which(grid.obj$logR > 2)
		grid.obj$logR[ii] <- jitter(rep(2, length(ii)),  amount=0.5)
	} else {
		ii <- which(grid.obj$logR > 1)
		if(length(ii) > 0){
			grid.obj$logR[ii] <- jitter(rep(1, length(ii)),  amount=0.5)
		}
	}
	plot.logr <- xyplot(logR ~ x | subject,
			    data=grid.obj,
			    pch=21,
			    cex=cex,
			    layout=c(1,4),
			    index.cond=list(4:1),
			    xlab="Mb",
			    ylab="",
			    ...)
	return(plot.logr)
}

grid.trioBaf <- function(mset, cex=0.3){
	grid.obj2 <- data.frame(baf=c(baf.F(mset),
				baf.M(mset),
				baf.O(mset)),
				x=rep(position(mset)/1e6, 3),
				subject=factor(rep(1:3, each=nrow(mset)), labels=c("Father", "Mother", "Offspring"), ordered=TRUE))
	plot.baf <- xyplot(baf ~ x | subject,
			   data=grid.obj2, pch=21, cex=cex,
			   layout=c(1,3), index.cond=list(3:1),
			   xlab="",
			   ylab="")
	return(plot.baf)
}

addSegments <- function(cbs.segs, distanceRanges, CHR, mset,...){
        ranges.md <- distanceRanges[distanceRanges$chrom == CHR, ]
	ranges.md$seg.mean <- -1*ranges.md$seg.mean
	if("ylim" %in% names(list(...))){
		ylim <- list(...)[["ylim"]]
		cbs.segs$seg.mean[cbs.segs$seg.mean < ylim[1]] <- ylim[1]+0.1
		ranges.md$seg.mean[ranges.md$seg.mean < ylim[1]] <- ylim[1]+0.1
	}
	trellis.focus("panel", 1, 4, highlight = FALSE)
	panel.grid(h = 10, v = 10, col = "grey", lty = 3)
	cbs.sub <- cbs.segs[cbs.segs$id==paste(sampleNames(mset)[1], "03", sep="_"), ]
	panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2)#gp=gpar("lwd"=2))

	trellis.focus("panel", 1, 3, highlight = FALSE)
	panel.grid(h = 10, v = 10, col = "grey", lty = 3)
	cbs.sub=cbs.segs[cbs.segs$id==paste(sampleNames(mset)[1], "02", sep="_"), ]
	panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2)#gp=gpar("lwd"=2))

	trellis.focus("panel", 1, 2, highlight = FALSE)
	panel.grid(h = 10, v = 10, col = "grey", lty = 3)
	cbs.sub <- cbs.segs[cbs.segs$id==paste(sampleNames(mset)[1], "01", sep="_"), ]
	tmp <- panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2)#gp=gpar("lwd"=2))

	cbs.sub <- ranges.md[substr(ranges.md$id, 1, 5) %in% sampleNames(mset), ]
	trellis.focus("panel", 1, 1, highlight = FALSE)
	panel.grid(h = 10, v = 10, col = "grey", lty = 3)
	tmp <- panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2)#gp=gpar("lwd"=2))
	upViewport(0)
}

addGenes <- function(rf, chr.name, xlim){
	rf.chr <- rf[rf$txStart <= xlim[2] & rf$txEnd >= xlim[1] & rf$chrom==chr.name, ]
	flatBed <- flatten.bed(rf.chr)
	panel.flatbed(flat=flatBed, showIso=FALSE, rows=5, cex=0.7)
}

addCnv <- function(cnv, chr.name, xlim, col="red", fill="red"){
	cnv.chr <- cnv[cnv$txStart <= xlim[2] & cnv$txEnd >= xlim[1] & cnv$chrom==chr.name, ]
	flatBed.cnv <- flatten.bed(cnv.chr)
	panel.flatbed(flat=flatBed.cnv, showIso=FALSE, rows=5, cex=0.7, col=col, fill=fill)
}

highlightRegion <- function(ranges.object, index,
			    col="grey80",
			    lwd=1,
			    bg="lightblue",
			    alpha=0.1){
	bg2 <- col2rgb(bg)/255
	fill <- rgb(bg2[1, ], bg2[2,], bg2[3,], alpha=alpha)
	for(k in seq_along(index)){
		ii <- index[k]
		st <- start(ranges.object[ii, ])/1e6
		en <- end(ranges.object[ii,])/1e6
		grid.rect(x=unit(st, "native"),
			  y=unit(0, "native"),
			  width=unit(en-st, "native"),
			  height=unit(1, "native"),
			  just=c("left", "bottom"), gp=gpar(col=col, lwd=lwd, fill=fill))
	}
}

trioIndices <- function(ranges.object, selectedRanges, xlim){
	CHR <- unique(selectedRanges$chrom)
	index <- which(ranges.object$chrom == CHR & ranges.object$id %in% selectedRanges$id & start(ranges.object) <= xlim[2] & end(ranges.object) >= xlim[1])
	index <- index[unlist(sapply(selectedRanges$id, function(x,y) grep(x, y), y=ranges.object$id[index]))]
	ids <- ranges.object$id[index]
	unique.ids <- unique(ids)
	index <- split(index, ids)
	index <- index[match(unique.ids, names(index))]
	index <- lapply(index, unique)
	return(index)
}

grid.trio <- function(ranges.object, distanceRanges, distanceSet,
		      bsSet,
		      index,
		      xlim,
		      CHR,
		      add.segments=TRUE,
		      highlight=TRUE,
		      alpha=0.1,
		      add.genes=TRUE,
		      cex.pch=0.3,
		      label="min distance", ...){
	chr.name <- paste("chr", CHR, sep="")
	for(j in seq_along(index)){
		grid.newpage()
		cat("Fig ", j, "\n")
		i <- index[[j]][1]
		mset <- constructTrioSetFromRanges(ranges.object[i, ],
						   distanceSet,
						   bsSet,
						   xlim=xlim)
		plot.logr <- grid.trioLogR(mset,
					   labels=c("Father", "Mother", "Offspring", label),
					   cex=cex.pch, ...)
		plot.baf <- grid.trioBaf(mset, cex=cex.pch)
		if(exists("cbs.segs")){
			if(CHR != unique(cbs.segs$chrom))
				load(file.path(beadstudiodir(), paste("cbs_chr" , CHR, ".rda", sep="")))
		} else 	load(file.path(beadstudiodir(), paste("cbs_chr" , CHR, ".rda", sep="")))
		st <- start(ranges.object[i, ])/1e6
		en <- end(ranges.object[i,])/1e6
		lvp <- viewport(x=0, width=unit(0.5, "npc"), just="left", name="lvp")
		pushViewport(lvp)
		pushViewport(dataViewport(xscale=plot.logr$x.limits, yscale=c(0,1), clip="on"))
		print(plot.logr, newpage=FALSE, prefix="plot1", more=TRUE)
		if(add.segments) addSegments(cbs.segs=cbs.segs, distanceRanges=distanceRanges, CHR=CHR, mset=mset, ...)
		if(highlight){
			pushViewport(viewport(x=unit(0.25, "npc"),
					      y=unit(0.5, "npc"),
					      width=unit(0.5, "npc"),
					      height=unit(1, "npc")))
			pushViewport(plotViewport(c(0, 3, 0, 3)))  ## needs to be exact same as above
			pushViewport(dataViewport(xscale=plot.logr$x.limits, yscale=c(0,1), clip="on"))
			suppressWarnings(highlightRegion(ranges.object, index[[j]], alpha=alpha))
		}
		upViewport(0)
		grid.text("Log R Ratio", x=unit(0.25, "npc"), y=unit(0.97, "npc"), gp=gpar("cex"=1.1))
		grid.text(paste(chr.name, ", Family", sampleNames(mset)), x=unit(0.5, "npc"), y=unit(0.98, "npc"), gp=gpar("cex"=1.2))

		## Plot BAFs
		upViewport(0)
		print(plot.baf, position=c(0.5, 0.21,1, 1), more=TRUE, prefix="baf")
		grid.text("B Allele Frequency", x=unit(0.75, "npc"), y=unit(0.97, "npc"), gp=gpar("cex"=1.1))
		if(add.genes){
			upViewport(0)
			pushViewport(viewport(x=unit(0.75, "npc"),
					      y=unit(0.12, "npc"),
					      width=unit(0.5, "npc"),
					      height=unit(0.12, "npc"), just="bottom"))
			pushViewport(plotViewport(c(0, 3, 0,3)))  ## needs to be exact same as above
			pushViewport(dataViewport(xscale=xlim/1000, yscale=c(0,1), clip="on"))
			addGenes(rf, chr.name, xlim)

			upViewport(0)
			pushViewport(viewport(x=unit(0.75, "npc"),
					      y=unit(0, "npc"),
					      width=unit(0.5, "npc"),
					      height=unit(0.12, "npc"), just="bottom"))
			pushViewport(plotViewport(c(0, 3, 0,3)))  ## needs to be exact same as above
			pushViewport(dataViewport(xscale=xlim/1000, yscale=c(0,1), clip="on"))
			addCnv(cnv, chr.name, xlim)
		}
		if(highlight){
			upViewport(0)
			pushViewport(viewport(x=unit(0.75, "npc"),
					      y=unit(0.5, "npc"),
					      width=unit(0.5, "npc"),
					      height=unit(1, "npc")))
			pushViewport(plotViewport(c(0, 3, 0,3)))  ## needs to be exact same as above
			pushViewport(dataViewport(xscale=plot.logr$x.limits, yscale=c(0,1), clip="on"))
			suppressWarnings(highlightRegion(ranges.object, index[[j]], alpha=alpha))
		}
	}
}

getXlimAmp <- function(chr){
	switch(chr,
	       chr1.all=c(0, 17e6),
	       chr1.1=c(0,5e6),
	       chr1.2=c(15e6, 17e6),
	       chr2=c(60e6, 62.5e6),
	       chr3.1=c(45e6, 55e6),
	       chr3.2=c(85e6, 98e6),
	       chr4=c(65e6, 75e6),
	       chr5=c(45e6, 52e6),
	       chr6=c(32.4e6, 33e6),
	       chr7=c(96e6, 105e6),
	       chr8=c(38e6, 39.8e6),
	       chr11=c(51.5e6, 58e6),
	       chr14=c(21.5e6, 22.5e6),
	       chr15=c(32.2e6, 33e6),
	       chr16=c(32e6, 33e6),
	       chr17=c(74.4e6, 75e6),
	       chr18=c(64e6, 66e6),
	       chr19=c(19.8e6, 21.4e6),
	       chr22=c(17.8e6, 19.8e6))
}


explainingMinDistance <- function(deletion.ranges,
				  show.segmentation=FALSE,
				  show.distance=FALSE,
				  show.distance.segmentation=FALSE,
				  highlight=FALSE,
				  label="min distance",
				  highlight.border="grey80",
				  highlight.col="lightblue",
				  highlight.alpha=0.1){
	xlim <- getGenomicAxis(deletion.ranges, CHR=CHR, FRAME=1e6)
	xlim[2] <- 22e6
	i <- grep("19225", names(index))
	i <- index[[i]][1]
	cex <- 0.3
	mset <- constructTrioSetFromRanges(deletion.ranges[i, ], minDistanceSet,
					   bsSet, xlim=xlim)
	grid.obj <- data.frame(logR=c(logR.F(mset),
			       logR.M(mset),
			       logR.O(mset),
			       -mindist(mset)),
			       x=rep(position(mset)/1e6, 4),
			       subject=factor(rep(1:4, each=nrow(mset)), labels=c("Father", "Mother", "Offspring", label), ordered=TRUE))
	ylim <- c(-1,0.7)
	plot.logr <- xyplot(logR ~ x | subject,
			    data=grid.obj,
			    pch=21,
			    cex=cex,
			    layout=c(1,4),
			    index.cond=list(4:1), ylim=ylim,
			    xlab="Mb",
			    ylab="",
			    panel=myPanel)
	lvp <- viewport(x=0, y=unit(1, "npc"),
			width=unit(1, "npc"), height=unit(1, "npc"),
			just=c("left", "top"), name="lvp")
	load(file.path(beadstudiodir(), paste("cbs_chr" , CHR, ".rda", sep="")))
	grid.newpage()
	pushViewport(lvp)
	print(plot.logr, newpage=FALSE, more=TRUE)
	if(!show.distance){
		trellis.focus("panel", 1, 1, highlight = FALSE)
		grid.points(grid.obj$x[grid.obj$subject==label],
			    grid.obj$logR[grid.obj$subject==label],
			    pch=21, gp=gpar(col="white", cex=0.5, fill="white"))
		panel.grid(h = 10, v = 10, col = "grey", lty = 3)
	}
	if(show.segmentation){
		trellis.focus("panel", 1, 4, highlight = FALSE)
		panel.grid(h = 10, v = 10, col = "grey", lty = 3)
		cbs.sub <- cbs.segs[cbs.segs$id==paste(sampleNames(mset)[1], "03", sep="_"), ]
		panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2)#gp=gpar("lwd"=2))
		trellis.focus("panel", 1, 3, highlight = FALSE)
		panel.grid(h = 10, v = 10, col = "grey", lty = 3)
		cbs.sub=cbs.segs[cbs.segs$id==paste(sampleNames(mset)[1], "02", sep="_"), ]
		panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2)#gp=gpar("lwd"=2))
		trellis.focus("panel", 1, 2, highlight = FALSE)
		panel.grid(h = 10, v = 10, col = "grey", lty = 3)
		cbs.sub <- cbs.segs[cbs.segs$id==paste(sampleNames(mset)[1], "01", sep="_"), ]
		panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2)#gp=gpar("lwd"=2))
	}
##	if(show.distance){
##		upViewport(0)
##		pushViewport(tvp)
##		plot(plot.logr2, position=c(newpage=FALSE)
##	}
	if(show.distance.segmentation){
		ranges.md <- minDistanceRanges[minDistanceRanges$chrom == CHR, ]
		ranges.md$seg.mean <- -1*ranges.md$seg.mean
		cbs.sub <- ranges.md[substr(ranges.md$id, 1, 5) %in% sampleNames(mset), ]
		trellis.focus("panel", 1, 1, highlight = FALSE)
		panel.grid(h = 10, v = 10, col = "grey", lty = 3)
		panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2)#gp=gpar("lwd"=2))
	}
	if(highlight){
		upViewport(0)
		pushViewport(viewport(x=unit(0.5, "npc"),
				      y=unit(0.5, "npc"),
				      width=unit(1, "npc"),
				      height=unit(0.98, "npc")))
		pushViewport(plotViewport(c(0, 3, 0, 3)))  ## needs to be exact same as above
		pushViewport(dataViewport(xscale=plot.logr$x.limits, yscale=c(0,1), clip="on"))
		i <- grep("19225", names(index))
		highlightRegion(deletion.ranges, index[[i]],
				col=highlight.border,
				bg=highlight.col,
				alpha=highlight.alpha)
	}
}

myboxplotPanel <- function(x,y, subscripts, col="black", cex=1, ...){
	panel.bwplot(x, y, col=col, cex=cex, ...)
	panel.grid(h = 10, v=0, col = "grey", lty = 3)
}
