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

