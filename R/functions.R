
projectMetadata <- function(path="/thumper/ctsa/beaty/holger/txtfiles", ## path to BeadStudioFiles
			    samplesheet,
			    outdir=beadstudiodir(),
			    bsSetFile="~/Projects/BeatyExperimentData/data/bsSetAll.rda"){
	fnames <- list.files(path, full.names=TRUE)


}

initializeBeadStudioContainer <- function(path="/thumper/ctsa/beaty/holger/txtfiles",
					  outdir=beadstudiodir(),
					  samplesheet,
					  cdfName="human610quadv1b",
					  filename="~/Projects/BeatyExperimentData/data/bsSetAll.rda",
					  verbose=TRUE){
	if(missing(samplesheet)) data(samplesheet)
	fnames <- list.files(path, full.names=TRUE)
	if(!file.exists(outdir)) dir.create(outdir)
	if(!file.exists(filename)){
		message(filename, " does not exist. Initializing bsSet")
		bsSet <- constructLogRatioSet(fnames, filename, cdfName, samplesheet)
		message("Saving bsSet...")
		save(bsSet, file=filename)
	}
	return(bsSet)
}
populateAssayData <- function(path, readScriptFilename="ReadBsScript2.R"){
	message("submitting batch jobs to populate assayData.  Must be on fat node")
	fnames <- list.files(path, full.names=TRUE)
	if(!file.exists(readScriptFilename)) stop(readScriptFilename, " not in working directory")
	ocSamples(500)
	batches <- splitIndicesByLength(seq(along=fnames), ocSamples())
	for(J in seq_along(batches)){
		sink("temp")
		cat("J <- ", J, "\n")
		cat("filename <- ", filename, "\n")
		sink()
		fn <- paste("rs_", J, ".R", sep="")
		if(file.exists(fn)) unlink(fn)
		system(paste("cat temp ReadBsScript2.R >", fn))
		system(paste('cat ~/bin/cluster.template | perl -pe "s/Rprog/rs_', J, '.R/" > rs_', J, '.R.sh', sep=""))
		system(paste("qsub -m e -r y -cwd -l mem_free=10G,h_vmem=16G rs_", J, ".R.sh", sep=""))
		Sys.sleep(60*3)
	}
}

readPedfile <- function(fnames, pedfile="~/Projects/Beaty/inst/extdata/may_peds.csv"){
	mayped <- read.csv(pedfile, sep=";", as.is=TRUE)
	is.father <- which(mayped$father != "0")
	is.mother <- which(mayped$mother != "0")
	## a complete trio would be the length of is.father and is.mother
	stopifnot(all.equal(is.father, is.mother))
	i <- is.father
	trio.matrix <- cbind(mayped$father[i],
			     mayped$mother[i],
			     mayped$cidr_name[i])
	rownames(trio.matrix) <- ss(trio.matrix[,1])
	index <- which(trio.matrix[, 1] %in% fnames & trio.matrix[, 2] %in% fnames & trio.matrix[, 3] %in% fnames)
	trio.matrix <- trio.matrix[index, ]
	colnames(trio.matrix) <- c("F", "M", "O")
	return(trio.matrix)
}

addPhenoData <- function(bsSet,
			 pedfile="~/Projects/Beaty/inst/extdata/may_peds.csv"){
	##
	## READ mayped file and merge with existing phenodata
	##
	##
	message("Adding additional annotation to phenoData")
	## add additional phenodata
	mayped <- read.csv(pedfile, sep=";", as.is=TRUE)
	is.father <- which(mayped$father != "0")
	is.mother <- which(mayped$mother != "0")
	## a complete trio would be the length of is.father and is.mother
	stopifnot(all.equal(is.father, is.mother))
	pd <- pData(bsSet)
	open(pd$MAD)
	pd2 <- merge(pd, mayped, by.x="CIDR_Name", by.y="cidr_name", all.x=TRUE)
	close(pd$MAD)
	stopifnot(all.equal(pd$CIDR_Name, pd2$CIDR_Name))
	tmp <- new("AnnotatedDataFrame", data=pd2)
	sampleNames(tmp) <- sampleNames(bsSet)
	phenoData(bsSet) <- tmp
	##
	## DNA source
	##
	dna <- getDnaSource(bsSet)
	wga.fam <- isWgaFamily(sssampleNames(bsSet), dna)
	bsSet$is.wga <- sssampleNames(bsSet)%in%wga.fam
	bsSet$dna <- getDnaSource(bsSet)
	##
	## who are the complete trios (exclude duplicates)
	##
	complete.trios <- completeTrios(bsSet)
	bsSet$complete.trio <- ssampleNames(bsSet) %in% as.character(complete.trios)
	##
	##
	##bsSet$complete.trio <- sssampleNames(bsSet) %in% complete.trios
	nTrios <- nrow(complete.trios)
	##nTrios <- length(complete.trios) ## 2065 F-M-O trios (duplicates removed)
	## Family name of trios flagged by QC
	flagged.trios <- qcFlag(bsSet)
	##
	##
	##
	## families not flagged by QC
	trio.names <- unique(ss(as.character(complete.trios)))
	trio.passqc <- trio.names[!trio.names %in% flagged.trios]
	##complete.trio.qc <- complete.trios[!complete.trios %in% flagged.trios] ##1694
	##
	##
	bsSet$pass.qc <- !(sssampleNames(bsSet) %in% flagged.trios)
	##
	## Duplicate Names
	##
	##
	library(xlsx)
	dat <- read.xlsx(file="/home/bst/student/rscharpf/Projects/Beaty/inst/extdata/duplicates_report.xlsx",1)
	duplicateNames <- as.character(dat$Sample.Name)
	is.duplicate <- sampleNames(bsSet) %in% paste(duplicateNames, ".txt", sep="")
	bsSet$is.duplicate <- is.duplicate
	##
	dat2 <- read.xlsx(file="/home/bst/student/rscharpf/Projects/Beaty/inst/extdata/duplicates_report.xlsx",2)
	dupName1 <- substr(as.character(dat2[[1]]), 18, 25)
	dupName2 <- substr(as.character(dat2[[2]]), 18, 25)
	##
	bsSet$duplicateSample <- NA
	index <- match(dupName1, ssampleNames(bsSet))
	bsSet$duplicateSample[index] <- dupName2
	##
	##
	index <- match(dupName2, ssampleNames(bsSet))
	bsSet$duplicateSample[index] <- dupName1
	##bsSet$duplicate.factor <- 0L
	##bsSet$duplicate.set <- 0L
	##ix1 <- match(dupName1, ssampleNames(bsSet))
	##ix2 <- match(dupName2, ssampleNames(bsSet))
	##bsSet$duplicate.factor[ix1] <- seq_along(dupName1)
	##bsSet$duplicate.set[ix1] <- 1L
	##bsSet$duplicate.factor[ix2] <- seq_along(dupName1)##+length(dupName1)
	##bsSet$duplicate.set[ix2] <- 2L
	##save(bsSet, file="~/Projects/BeatyExperimentData/data/bsSetAll.rda")
	##save(bsSet, file=filename)
	return(bsSet)
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

readProcessedFiles <- function(bsSet, stratum){
	sns <- sampleNames(bsSet)
	path <- "/thumper/ctsa/beaty/holger/txtfiles"
	fnames <- list.files(path, full.names=TRUE)
	names(fnames) <- substr(basename(fnames), 1, 8)
	table(names(fnames) %in% sns)
	fnames <- fnames[names(fnames) %in% sns]
	fnames <- fnames[match(sns, names(fnames))]
	stopifnot(identical(names(fnames), sns))

	dat <- read.delim(fnames[1], colClasses=c("character", "numeric", "numeric"))
	match.index <- match(featureNames(bsSet), dat$Name)
	j <- splitIndicesByLength(seq(along=fnames), ocSamples())[[stratum]]

	tmpBaf <- matrix(NA, nrow(bsSet), length(j))
	tmpLogr <- matrix(NA, nrow(bsSet), length(j))
	for(k in seq(along=j)){
		if(k %% 10 == 0) cat(".")
		if(k > 1){
			dat <- read.delim(fnames[j[k]], colClasses=c("character", "numeric", "numeric"))
		}
		dat <- dat[match.index, ]
		if(k == 1) stopifnot(all.equal(dat$Name, featureNames(bsSet)))
		##order baf and logr values by chromosome, physical position
		tmpBaf[, k] <- dat[, 3]
		tmpLogr[, k] <- dat[, 2]
	}
	autosome.index <- which(chromosome(bsSet) < 23)
	mads <- apply(tmpLogr, 2, mad, na.rm=TRUE)
	open(baf(bsSet))
	open(logR(bsSet))
	baf(bsSet)[, j] <- tmpBaf
	logR(bsSet)[, j] <- tmpLogr
	bsSet$MAD[j] <- as.numeric(mads)
	close(baf(bsSet))
	close(logR(bsSet))
	TRUE
}

readProcessedFiles2 <- function(bsSet, stratum){
	sns <- sampleNames(bsSet)
	path <- "/thumper/ctsa/beaty/holger/txtfiles"
	fnames <- list.files(path, full.names=TRUE)
	##names(fnames) <- substr(basename(fnames), 1, 8)
	names(fnames) <- basename(fnames)
	##table(names(fnames) %in% sns)
	##fnames <- fnames[names(fnames) %in% sns]
	fnames <- fnames[match(sns, names(fnames))]
	stopifnot(identical(names(fnames), sns))

	dat <- read.delim(fnames[1], colClasses=c("character", "numeric", "numeric"))
	match.index <- match(featureNames(bsSet), dat$Name)
	j <- splitIndicesByLength(seq(along=fnames), ocSamples())[[stratum]]

	tmpBaf <- matrix(NA, nrow(bsSet), length(j))
	tmpLogr <- matrix(NA, nrow(bsSet), length(j))
	for(k in seq(along=j)){
		if(k %% 10 == 0) cat(".")
		if(k > 1){
			dat <- read.delim(fnames[j[k]], colClasses=c("character", "numeric", "numeric"))
		}
		dat <- dat[match.index, ]
		if(k == 1) stopifnot(all.equal(dat$Name, featureNames(bsSet)))
		##order baf and logr values by chromosome, physical position
		tmpBaf[, k] <- dat[, 3]
		tmpLogr[, k] <- dat[, 2]
	}
	autosome.index <- which(chromosome(bsSet) < 23)
	mads <- apply(tmpLogr, 2, mad, na.rm=TRUE)
	open(baf(bsSet))
	open(logR(bsSet))
	baf(bsSet)[, j] <- tmpBaf
	logR(bsSet)[, j] <- tmpLogr
	bsSet$MAD[j] <- as.numeric(mads)
	close(baf(bsSet))
	close(logR(bsSet))
	TRUE
}

cbsSegmentation <- function(object, sample.index, alpha.parent=0.05, alpha.offspring=0.01){
	## allow segmentation on parents to jump more easily than for offspring
	parent.index <- sample.index[object$pedId[sample.index]=="father" | object$pedId[sample.index]=="mother"]
	offspring.index <- sample.index[object$pedId[sample.index]=="offspring"]
	marker.index <- which(chromosome(object) == CHR & !duplicated(position(object)))
	open(logR(object))
	CNA.object <- CNA(genomdat=as.matrix(logR(object)[marker.index, parent.index]),
			  chrom=chromosome(object)[marker.index],
			  maploc=position(object)[marker.index],
			  data.type="logratio",
			  sampleid=sampleNames(object)[parent.index])
	smu.object <- smooth.CNA(CNA.object)
	tmp <- segment(smu.object, verbose=0, alpha=0.05)
	cbs.segs1 <- cbind(tmp$output, tmp$segRows)
	CNA.object <- CNA(genomdat=as.matrix(logR(object)[marker.index, offspring.index]),
			  chrom=chromosome(object)[marker.index],
			  maploc=position(object)[marker.index],
			  data.type="logratio",
			  sampleid=sampleNames(object)[offspring.index])
	smu.object <- smooth.CNA(CNA.object)
	tmp <- segment(smu.object, verbose=0, alpha=0.01)
	cbs.segs2 <- cbind(tmp$output, tmp$segRows)
	cbs.segs <- rbind(cbs.segs1, cbs.segs2)
	cbs.segs
}

constructMinDistanceContainer <- function(bsSet){
	trios <- completeTrios(bsSet)
	##offspring.index <- which(who(sampleNames(bsSet))=="offspring")
	offspring.index <- match(trios[, "O"], ssampleNames(bsSet))
	ldPath(beadstudiodir())
	min.dist <- initializeBigMatrix("min.dist2", nr=nrow(bsSet), nc=length(offspring.index), vmode="double")
	##cnConf <- initializeBigMatrix("cnConf2", nr=nrow(bsSet), nc=length(offspring.index), vmode="integer")
	dimnames(min.dist) <- list(featureNames(bsSet), sampleNames(bsSet)[offspring.index])
	minDistanceSet <- new("MinDistanceSet",
			      mindist=min.dist,
			      ##cnConfidence=cnConf,
			      phenoData=phenoData(bsSet)[offspring.index, ],
			      annotation=annotation(bsSet))
	mads <- initializeBigVector(name="mindist.mads", n=ncol(minDistanceSet), vmode="double")
	minDistanceSet$MAD <- mads
	featureData(minDistanceSet) <- featureData(bsSet)
	minDistanceSet
}

segmentMD <- function(minDistanceSet,
		      id,
		      verbose=FALSE, ...){
	## needs to be ordered
	ix <- order(chromosome(minDistanceSet), position(minDistanceSet))
	stopifnot(all(diff(ix) > 0))
	##
	##
	dfl <- vector("list", 22) ## data.frame list
	ix <- match(s(id), ssampleNames(minDistanceSet))
	stopifnot(length(ix) > 0)
	open(minDistanceSet)
	##
	##
	marker.index.list <- split(seq(length=nrow(minDistanceSet)), chromosome(minDistanceSet))
	for(CHR in 1:22){
		dfl[[CHR]] <- segmentBatchWithCbs(minDistanceSet=minDistanceSet,
						  marker.index=marker.index.list[[CHR]],
						  sample.index=ix,
						  CHR=CHR,
						  verbose=verbose)
	}
	close(minDistanceSet)
	df <- do.call("rbind", dfl)
	return(df)
}

segmentBatchWithCbs <- function(minDistanceSet,
				marker.index,
				sample.index,
				CHR,
				verbose=FALSE, ...){
	if("undo.splits" %in% names(list(...))) message("undo.splits = ", list(...)[["undo.splits"]])
	stopifnot(!missing(CHR))
	fns <- featureNames(minDistanceSet)  ## Do not subset!!
	##
	pos <- position(minDistanceSet)[marker.index]
	marker.index <- marker.index[!duplicated(pos)]
	pos <- position(minDistanceSet)[marker.index]
	chrom <- chromosome(minDistanceSet)[marker.index]
	CN <- copyNumber(minDistanceSet)[marker.index, sample.index, drop=FALSE]
	id <- sampleNames(minDistanceSet)[sample.index]
	arm <- getChromosomeArm(chrom, pos)
	index.list <- split(seq_along(marker.index), arm)
	md.segs <- list()
	if(verbose) message("Running CBS by chromosome arm")
	for(i in seq_along(index.list)){
		j <- index.list[[i]]
		CNA.object <- CNA(genomdat=CN[j, , drop=FALSE],
				  chrom=chrom[j],
				  maploc=pos[j],
				  data.type="logratio",
				  sampleid=id)
		smu.object <- smooth.CNA(CNA.object)
		tmp <- segment(smu.object, verbose=0, ...)
		df <- tmp$output
		sr <- tmp$segRows
		##df <- cbind(tmp$output, tmp$segRows)
		##md.segs[[i]] <-
		firstMarker <- rownames(CNA.object)[sr$startRow]
		endMarker <- rownames(CNA.object)[sr$endRow]
		df$start.index <- match(firstMarker, fns)
		df$end.index <- match(endMarker, fns)
		md.segs[[i]] <- df
	}
	if(length(md.segs) > 1){
		md.segs <- do.call("rbind", md.segs)
	} else md.segs=md.segs[[1]]
	## md.segs is a data.frame
	md.segs
}

submitCbsJobs <- function(BATCHSIZE, NN, CHR, undo.splits="'none'", outdir){
	for(BATCH in 1:NN){
		sink("tempMD")
		cat("BATCH <- ", BATCH, "\n")
		cat("BATCHSIZE <- ", BATCHSIZE, "\n")
		cat("CHR <- ", CHR, "\n")
		##cat("undo.splits <-", undo.splits, "\n")
		cat("outdir <- '", outdir, "'\n", sep="")
		sink()
		fn <- paste("mds_c", CHR, "_", BATCH, ".R", sep="")
		if(file.exists(fn)) unlink(fn)
		system(paste("cat tempMD RunMinDistanceSegmentation.R >", fn))
		system(paste('cat ~/bin/cluster.template | perl -pe "s/Rprog/mds_c', CHR, "_", BATCH, '.R/" > mds_c', CHR, "_", BATCH, '.R.sh', sep=""))
		system(paste("qsub -m e -r y -cwd -l mem_free=2G,h_vmem=3G mds_c", CHR, "_", BATCH, ".R.sh", sep=""))
		Sys.sleep(60*1) ## keeps multiple jobs from trying to load minDistanceSet simultaneously
	}
	return(NULL)
}

wrapperPosteriorProbs <- function(filename="~/Projects/Beaty/inst/scripts/BayesFactor.Rnw",
				  chromosomes, G=10, sleep=FALSE){
	stopifnot(file.exists(filename))
	Stangle(filename)
	for(i in seq_along(chromosomes)){
		CHR <- chromosomes[i]
		sink("temp")
		cat("CHR <- ", CHR, "\n")
		sink()
		fn <- paste("bf_", CHR, ".R", sep="")
		if(file.exists(fn)) unlink(fn)
		system(paste("cat temp BayesFactor.R >", fn))
		system(paste('cat ~/bin/cluster.template | perl -pe "s/Rprog/bf_', CHR, '.R/" > bf_', CHR, '.R.sh', sep=""))
		system(paste("qsub -m e -r y -cwd -l mem_free=", G, "G,h_vmem=",G+3,"G bf_", CHR, ".R.sh", sep=""))
		##if(CHR %% 15 == 0) Sys.sleep(60*60*5) ## keeps multiple jobs from trying to load minDistanceSet simultaneously
	}
	if(sleep) Sys.sleep(60*60*5)##3 minutes
	return(TRUE)
}

wrapperPrune <- function(chromosomes){
	Stangle("~/Projects/Beaty/inst/scripts/Prune.Rnw")
	for(i in seq_along(chromosomes)){
		CHR <- chromosomes[i]
		sink("temp")
		cat("CHR <- ", CHR, "\n")
		sink()
		fn <- paste("Prune_", CHR, ".R", sep="")
		if(file.exists(fn)) unlink(fn)
		system(paste("cat temp Prune.R >", fn))
		system(paste('cat ~/bin/cluster.template | perl -pe "s/Rprog/Prune_', CHR, '.R/" > Prune_', CHR, '.R.sh', sep=""))
		system(paste("qsub -m e -r y -cwd -l mem_free=10G,h_vmem=16G Prune_", CHR, ".R.sh", sep=""))
		##if(CHR %% 15 == 0) Sys.sleep(60*60*5) ## keeps multiple jobs from trying to load minDistanceSet simultaneously
	}
	##Sys.sleep(60*60*5)##3 minutes
	return(TRUE)
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

loadRangesCbs <- function(outdir, pattern, CHR, name){
	fname <- list.files(outdir, pattern=pattern, full.name=TRUE)
	if(missing(name)) stop("must specify object name")
	stopifnot(length(fname) == 1)
	load(fname)
	cbs.segs <- get(name)
	cbs.segs
}


##loadRanges <- function(outdir, pattern, CHR, name){
readCbsBatchFiles <- function(outdir, pattern, CHR, name){
	fnames <- list.files(outdir, pattern=pattern, full.name=TRUE)
	if(missing(name)) stop("must specify R object name when saved.")
	if(length(fnames) == 0) stop(paste("There are no segmentation files for chrom", CHR))
	segmeans <- vector("list", length(fnames))
	for(i in seq_along(segmeans)){
		load(fnames[i])
		cbs.segs <- get(name)
##		rd <- RangedData(IRanges(start=cbs.segs$startRow,
##					 end=cbs.segs$endRow),
##				 pos.start=cbs.segs$loc.start,
		rd <- RangedData(IRanges(cbs.segs$loc.start,
					 cbs.segs$loc.end),
				 ##pos.end=cbs.segs$loc.end,
				 ##num.mark=cbs.segs$num.mark,
				 id=cbs.segs$ID,
				 num.mark=cbs.segs$num.mark,
				 chrom=cbs.segs$chrom,
				 seg.mean=cbs.segs$seg.mean)
		segmeans[[i]] <- rd
	}
	rdlist <- RangedDataList(segmeans)
	rd <- stack(rdlist)
	ix <- match("sample", colnames(rd))
	if(length(ix) > 0) rd <- rd[, -ix]
	rd$id <- substr(rd$id, 2, 9)
	return(rd)
##	segmean_ranges <- do.call("c", segmeans)
##	if(substr(segmean_ranges$id[1], 1, 1) == "X"){
##		segmean_ranges$id <- substr(segmean_ranges$id, 2, 9)
##	}
##	segmean_ranges <- RangedData(IRanges(segmean_ranges$pos.start,
##					     segmean_ranges$pos.end),
##				     space=rep(CHR, nrow(segmean_ranges)),
##				     id=segmean_ranges$id,
##				     chrom=segmean_ranges$chrom,
##				     num.mark=width(segmean_ranges), ## would be number of markers that were not NAs (I think)
##				     seg.mean=segmean_ranges$seg.mean,
##				     pedId=who(segmean_ranges$id))
##	segmean_ranges
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

featuresInRange <- function(object, range, ...){
	##if(is(range, "RangedData")) stopifnot(nrow(range)==1) else  stopifnot(length(range)==1)
	##if(is(range, "GRanges")) CHR <- chromosome(range) else CHR <- range$chrom
	##featuresInXlim(object, start=start(range), end=end(range), CHR=CHR, ...)
}

addIndicesFromFeatureData <- function(rd.object, fD, FRAME=0){
	uchrom.rd <- unique(rd.object$chrom)
	uchrom.fd <- unique(fD$chromosome)
	stopifnot(length(uchrom.rd)==1)
	stopifnot(length(uchrom.fd)==1)
	stopifnot(uchrom.rd == uchrom.fd)
	##rd.object$start.index <- NA
	##rd.object$end.index <- NA

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

cbsLoader <- function(segname, CHR, envir){
	theFile <- file.path(beadstudiodir(), paste("cbs_chr" , CHR, ".rda", sep=""))
	beatyPkgLoader(segname, theFile, "cbs.segs", envir)
}
beatyPkgLoader <- function(objectname, filename, loadedName, envir){
	needToLoad <- !crlmm:::isLoaded(objectname, envir)
	if(needToLoad) {
		load(filename)
		.object <- get(loadedName)
		message("assigning ", objectname, " to environment", environmentName(envir))
		assign(objectname, .object, envir=envir)
	} else  message("cbs.segs object already loaded")
	crlmm:::getVarInEnv(objectname, environ=envir)
}


readPennCnv <- function(penndir="/thumper/ctsa/beaty/holger/penncnv/jointDat"){
	if(length(grep("trioDat", penndir))==1) warning("this is not the recommended aproach")
	fnames <- list.files(penndir, full.names=TRUE)
	tmp <- read.delim(fnames[1], nrows=5, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	colClasses <- c("integer", "integer", "integer", "integer", "character", "character",
			"character", "character", "character",
			"character")
	rdList <- vector("list", length(fnames))
	message("Reading ", length(rdList), " files")
	for(i in seq_along(fnames)){
		if(i %% 100 == 0) cat('.')
		penn.joint <- read.delim(fnames[i], header=TRUE, sep="\t", stringsAsFactors=FALSE,
					 colClasses=colClasses)
		penn.joint$LengthCNV <- as.integer(gsub(",", "", penn.joint$LengthCNV))
		chr <- paste("chr", penn.joint$Chromosome, sep="")
		rd <- RangedData(IRanges(penn.joint$StartPosition,
					 penn.joint$EndPosition),
				 chrom=penn.joint$Chromosome,
				 num.mark=penn.joint$NumberSNPs,
				 id=penn.joint$ID,
				 triostate=penn.joint$TrioState)
				 ##pedId=penn.joint$FamilyMember,
				 ##filename=basename(fnames[i]))
		rdList[[i]] <- rd
	}
	rdl <- RangedDataList(rdList)
	rd <- stack(rdl)
	ix <- match("sample", colnames(rd))
	if(length(ix) > 0) rd <- rd[, -ix]
##	penn.joint <- do.call("c", rdList)
##	if(offspring.only){
##		message("Returning only the ranges for the offspring")
##		##penn.joint <- penn.joint[penn.joint$pedId=="offspring", ]
##		fmo <- substr(penn.joint$id, 8, 8)
##		penn.joint <- penn.joint[fmo == 1, ]
##		##chrom <- as.character(space(penn.joint))
##		##chrom <- substr(chrom, 4, nchar(chrom))
##		##penn.joint$chrom <- chromosome2integer(chrom)
##	}
##	penn.joint2 <- RangedData(IRanges(start(penn.joint), end(penn.joint)),
##				  chrom=penn.joint$chrom,
##				  num.mark=penn.joint$num.mark,
##				  id=penn.joint$id,
##				  triostate=penn.joint$triostate,
##				  pedId=penn.joint$pedId,
##				  filename=penn.joint$filename)
##	return(penn.joint2)
	return(rd)
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

askHolger <- function(trios){
	trio.ids <- as.character(trios)
	I <- !(trio.ids %in%  sampleNames(penn.all))
	I <- matrix(I, ncol=3)
	index <- which(rowSums(I) > 0)
	incomplete <- trios[index, ]
	if(FALSE) save(incomplete, file="~/Projects/Beaty/data/askHolger.rda")
	return(incomplete)
}

pennStats <- function(penn.all, penn.offspring){
	id <- penn.all$id
	nSamples.penn <- length(unique(id))
	nTrios.penn <- length(unique(ss(id)))
	stopifnot(nSamples.penn/3 == nTrios.penn)
	triostates <- penn.all$triostate
	MIN.COV <- 10
	##fmo <- substr(penn.all$id, 8, 8)
	nM <- nMarkers(penn.all)
	is.o <- sampleNames(penn.all) %in% trios[, "O"]

	del.o <- which(substr(triostates, 3, 3) < 3 & nM >=MIN.COV & is.o)
	amp.o <- which(substr(triostates, 3, 3) > 3 & nM >=MIN.COV & is.o)
	ndel.all <- sum(substr(triostates, 3, 3) < 3 & is.o)
	namp.all <- sum(substr(triostates, 3, 3) > 3 & is.o)
	ndel.o <- length(del.o)
	namp.o <- length(amp.o)
	avgndel.o <- median(sapply(split(del.o, id[del.o]), length))
	avgnamp.o <- median(sapply(split(amp.o, id[amp.o]), length))
	cov.o <- median(nM[del.o])
	## for offspring
	nMo <- nMarkers(penn.offspring)
	tstate <- penn.offspring$triostate
	id <- sampleNames(penn.offspring)
	denovo.hemi <- which(tstate %in% Beaty:::offspring.hemizygous() & nMo >= MIN.COV)
	ndenovo.hemi <- length(denovo.hemi)
	avgndenovo.hemi <- median(sapply(split(denovo.hemi, id[denovo.hemi]), length))
	denovo.hom <- which(tstate %in% Beaty:::offspring.homozygous() & nMo >= MIN.COV)
	ndenovo.hom <- length(denovo.hom)
	avgndenovo.hom <- median(sapply(split(denovo.hom, id[denovo.hom]), length))
	denovo.dup <- which(tstate %in% Beaty:::duplicationStates() & nMo >= MIN.COV)
	ndenovo.dup <- length(denovo.dup)
	avgndenovo.dup <- median(sapply(split(denovo.dup, id[denovo.dup]), length))
	stats <- c(ndel.all,
		   namp.o,
		   ndel.o,
		   namp.all,
		   avgndel.o,
		   avgnamp.o,
		   ndenovo.hemi,
		   ndenovo.dup,
		   ndenovo.hom)
	names(stats) <- c("ndel.all",
			  "namp.o",
			  "ndel.o",
			  "namp.all",
			  "avgndel.o",
			  "avgnamp.o",
			  "ndenovo.hemi",
			  "ndenovo.dup",
			  "ndenovo.hom")
	return(stats)
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

collectAllRangesOfSize <- function(SIZE, bsSet,
				   minDistanceSet,
				   minDistanceRanges,
				   outdir, MIN=1, MAX=4, lambda=0.1,
				   MIN.THR=0.1,
				   denovo.deletion=TRUE){
	if(denovo.deletion){
		ranges2 <- minDistanceRanges[minDistanceRanges$seg.mean < 0 & minDistanceRanges$num.mark >= SIZE, ]
	} else {
		ranges2 <- minDistanceRanges[minDistanceRanges$seg.mean > 0 & minDistanceRanges$num.mark >= SIZE, ]
	}
	index <- match(ranges2$id, sampleNames(minDistanceSet))
	open(minDistanceSet$MAD)
	mads <- minDistanceSet$MAD[index]
	coverage <- ranges2$num.mark
	p <- lambda*exp(-lambda*coverage)
	b <- 1/(MAX - MIN)
	a <- MIN * b
	numberMads <- ((p-min(p))/(max(p)-min(p)) + a)/b
	thr <- numberMads * mads
	thr[thr < MIN.THR] <- MIN.THR
	if(denovo.deletion){
		thr <- -1*thr
		is.altered <- ifelse(ranges2$seg.mean <= thr, TRUE, FALSE)
		altered.RD <- ranges2[is.altered, ]
	} else {
		is.altered <- ifelse(ranges2$seg.mean >= thr, TRUE, FALSE)
		altered.RD <- ranges2[is.altered, ]
	}
	altered.RD
}

getFMOindex <- function(bsSet){
	family.id <- unique(substr(sampleNames(bsSet), 1, 5))
	FMO.index <- matrix(NA, length(family.id), 3)
	dimnames(FMO.index) <- list(family.id, c("F", "M", "O"))
	F.index <- match(paste(family.id, "03", sep="_"), substr(sampleNames(bsSet), 1, 8))
	M.index <- match(paste(family.id, "02", sep="_"), substr(sampleNames(bsSet), 1, 8))
	O.index <- match(paste(family.id, "01", sep="_"), substr(sampleNames(bsSet), 1, 8))
	FMO.index[, "F"] <- F.index
	FMO.index[, "M"] <- M.index
	FMO.index[, "O"] <- O.index
	FMO.index
}

rangeStats <- function(ranges.object, bsSet){
	pHet.2 <- pHet <- rep(NA, nrow(ranges.object))
	FMO.median <- matrix(NA, nrow(ranges.object), 3)
	colnames(FMO.median) <- c("F", "M", "O")
	F.index <- match(paste(substr(ranges.object$id, 1, 5), "03", sep="_"), substr(sampleNames(bsSet), 1, 8))
	M.index <- match(paste(substr(ranges.object$id, 1, 5), "02", sep="_"), substr(sampleNames(bsSet), 1, 8))
	O.index <- match(paste(substr(ranges.object$id, 1, 5), "01", sep="_"), substr(sampleNames(bsSet), 1, 8))
	FMO.index <- cbind(F.index, M.index, O.index)
	##sample.index <- match(ranges.object$id, sampleNames(bsSet))
	uid <- unique(ranges.object$id)
	for(j in seq_along(uid)){
		if(j %% 100 == 0) cat(".")
		id <- uid[j]
		fmo.i <- FMO.index[match(id, ranges.object$id), ]
		b <- as.numeric(baf(bsSet)[, fmo.i[3]])
		L <- as.matrix(logR(bsSet)[, fmo.i])
		range.index <- which(ranges.object$id == id)

		I <- cbind(ranges.object$start.index[range.index],
			   ranges.object$end.index[range.index])
		interp <- apply(I, 1, function(x) ":"(x[1], x[2]))
		fact <- rep(1:length(interp), sapply(interp, length))

		## proportion heterozygous
		i <- unlist(interp)
		bb <- b[i]
		bb.sp <- split(bb, fact)
		p <- sapply(bb.sp, function(x) mean(x >= 0.3 & x <= 0.7, na.rm=TRUE))
		p2 <- sapply(bb.sp, function(x) mean(x >= 0.45 & x <= 0.55, na.rm=TRUE))
		pHet[range.index] <- p
		pHet.2[range.index] <- p2
		## median log R Ratio
		LL <- L[i, ]
		for(jj in 1:3){
			ll.sp <- split(LL[, jj], fact)
			mn <- sapply(ll.sp, median, na.rm=TRUE)
			FMO.median[range.index, jj] <- mn
		}
	}
	return(list(pHet=pHet, pHet.2=pHet.2, FMO.median=FMO.median))
}

##			for(i in seq_along(index.list)){
##				ii <- match(fns.list[[i]], rownames(B))
##				b <- B[ii, i]
##				pHet[i] <- mean(b > 0.45 & b < 0.55, na.rm=TRUE)
##			}

## for denovo-amplifications
##collectAllRangesOfSize2 <- function(SIZE, bsSet,
##				    maxDistanceSet,
##				    maxDistanceRanges,
##				    outdir, MIN=1, MAX=4,
##				    lambda=0.1, upper.limit=-0.5){
##	data(chromosomeAnnotation)
##	centromere.ranges <- GRanges(seqnames=Rle(paste("chr", 1:22, sep=""), rep(1,22)),
##				     ranges=IRanges(chromosomeAnnotation[1:22, "centromereStart"],
##				     chromosomeAnnotation[1:22, "centromereEnd"]))
##	##ranges2 <- maxDistanceRanges[maxDistanceRanges$seg.mean < 0 & maxDistanceRanges$num.mark >= SIZE, ]
##	ranges2 <- maxDistanceRanges[maxDistanceRanges$seg.mean < upper.limit & maxDistanceRanges$num.mark >= SIZE, ]
##	index <- match(ranges2$id, sampleNames(maxDistanceSet))
##	open(maxDistanceSet$MAD)
##	mads <- maxDistanceSet$MAD[index]
##	x <- ranges2$num.mark
##	p <- lambda*exp(-lambda*x)
##	MIN <- 1; MAX <- 4
##	b <- 1/(MAX - MIN)
##	a <- MIN * b
##	numberMads <- ((p-min(p))/(max(p)-min(p)) + a)/b
##	thr <- -numberMads * mads
##	thr[thr > upper.limit] <- upper.limit
##	ranges2$is.altered <- ifelse(ranges2$seg.mean <= thr, TRUE, FALSE)
##	altered.RD <- ranges2[ranges2$is.altered, ]
##	altered.RD <- altered.RD[order(altered.RD$chrom, start(altered.RD)), ]
##	return(altered.RD)
##}

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


pruneByFactor <- function(range.object, f){
	rd <- list()
	for(i in seq_along(unique(range.object$id))){
		if(i %% 100==0) message(i, "/", length(unique(range.object$id)))
		id <- unique(range.object$id)[i]
		(index <- which(range.object$id == id))
		##trace(combineRangesByFactor, browser)
		rd[[i]] <- combineRangesByFactor(range.object[index, ], f=f[index])
	}
	ok <- tryCatch(tmp <- do.call("rbind", rd), error=function(e) FALSE)
	if(!ok) tmp <- rd
	return(tmp)
}

combineRangesByFactor <- function(range.object, f){
	ff <- cumsum(c(0, abs(diff(as.integer(as.factor(f))))))
	if(!any(duplicated(ff))) return(range.object)
	for(i in seq_along(unique(ff))){
		x <- unique(ff)[i]
		if(sum(ff==x) == 1) next()
		index <- which(ff==x)
		min.index <- min(index)
		max.index <- max(index)
		end(range.object)[index] <- max(end(range.object)[index])
		range.object$bayes.factor[index] <- sum(range.object$bayes.factor[index])
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
##
##
##		if(length(j) == 0){
##			j <- 1
##		}  else {
##		}
##
##		j <- j[-index]
##
##
##		j <- j[-index]
##
##
##		if(min.index > 1){
##			index.before <- 1:(min.index-1)
##		}
##		if(nrow(range.object) > max.index){
##			index.after <- (max.index+1):nrow(range.object)
##		}
##		if(min.index > 1 & nrow(range.object) > max.index){
##			index.new <- c(index.before, min.index, index.after)
##		}
##		if(min.index > 1 & nrow(range.object) <= max.index){
##			index.new <- c(index.before, min.index)
##		}
##		if(min.index == 1 & nrow(range.object) > max.index){
##			index.new <- c(min.index, index.after)
##		}
##		if(min.index == 1 & nrow(range.object) <= max.index){
##			index.new <- min.index
##		}
		range.object <- range.object[j, ]
	}
	return(range.object)
}

combine.data.frames <- function(dist.df, penn.df){
	if(is.null(dist.df) & is.null(penn.df)) return(NULL)
	if(is.null(dist.df)) dist.df <- penn.df[integer(0), ]
	if(is.null(penn.df)) penn.df <- dist.df[integer(0), ]
	##dist.df$method <- rep("distance", nrow(dist.df))
	##penn.df$method <- rep("penncnv", nrow(penn.df))
	combined.df <- rbind(dist.df, penn.df)
	combined.df <- combined.df[order(combined.df$chr), ]
##	maxim.chr <- split(combined.df$x1, combined.df$chr)
##	maxim <- sapply(maxim.chr, max)
##	maxim <- rep(maxim, sapply(maxim.chr, length))
##	min.chr <- split(combined.df$x0, combined.df$chr)
##	minim <- sapply(min.chr, min)
##	minim <- rep(minim, sapply(min.chr, length))
##	combined.df$min <- minim
##	combined.df$max <- maxim
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



harmonizeStates <- function(penn.joint, filter.multistate=FALSE){
	## note: 221 is denovo in the sense that neither parent had a homozygous deletion
	del.states <- deletionStates()
	amp.states <- duplicationStatesPenn()
	alt.states <- c(del.states, amp.states)
	##if(!all(penn.joint$pedId == "offspring")) stop("only offspring ranges can be in the object")
	##stopifnot(all(substr(penn.joint$id, 8, 8) == "1"))
	##colnames(penn.joint)[2] <- "num.mark"
	##penn.joint <- penn.joint[, -5] ## redundant
	message("Treating LOH state as 'normal'")
	penn.joint$state <- penn.joint$triostate
	penn.joint$state <- gsub("4", "3", penn.joint$state)
	message("Substituting '4' for states 5 and 6'")
	penn.joint$state <- gsub("5", "4", penn.joint$state)
	penn.joint$state <- gsub("6", "4", penn.joint$state)
	index.multiple.states <- grep("-", penn.joint$state)
	if(!filter.multistate){
		multi.state <- penn.joint$state[index.multiple.states]
		if(length(multi.state) > 0){
			message("Several regions have multiple states assigned... returning the first denovo state in the list")
			##multi.state <- sapply(multi.state, function(x) unique(strsplit(x, "-")[[1]]))
			checkMultiState <- function(x){
				state <- strsplit(x, "-")[[1]]
				if(any(state %in% alt.states)){
					state <- state[which(state %in% alt.states)[1]]
				} else{
					state <- state[1]
				}
				return(state)
			}
			state.cat <- sapply(multi.state, checkMultiState)
			##penn.joint$triostate[index.multiple.states] <- state.cat
			penn.joint$state[index.multiple.states] <- state.cat
		}
	} else message("multi-state ranges left as is")
	##penn.joint <- penn.joint[penn.joint$triostate %in% alt.states, ] ##68,472
	##return(penn.joint)
	penn.joint$state
}



callRange <- function(rd.object, bsSet){
	if(!"pedId" %in% varLabels(bsSet)) bsSet$pedId <- who(sampleNames(bsSet))
	##stopifnot(nrow(rd.object) == 1)
	##CHR <- rd.object$chrom[[1]]
	##marker.index <- which(chromosome(bsSet) == CHR & position(bsSet) >= start(rd.object) & position(bsSet) <= end(rd.object))
	family <- substr(sampleNames(bsSet), 1, 5)
	## each range corresponds to one trio
	##
	##sample.index <- which(family == substr(rd.object$id, 1, 5))
	##pedId <- bsSet$pedId[sample.index]
	##sample.index <- sample.index[match(pedId, c("father", "mother", "offspring"))]
	meds <- apply(logR(bsSet)[marker.index, sample.index], 2, "median", na.rm=T)
	homo.deletion <- meds <= -1
	hemi.deletion <- meds < -0.2 & meds > -1
	normal <- meds > -0.2 & meds < 0.1
	amp <- meds > 0.1
	calls <- rbind(homo.deletion, hemi.deletion, normal, amp)
	state <- paste(apply(calls, 2, which), collapse="")
	return(state)
}

##




calculateTrioState <- function(ranges.object, bsSet){
	state <- rep(NA, nrow(ranges.object))
	for(i in 1:nrow(ranges.object)){
		if(i %% 100 == 0) cat(".")
		state[i] <- callRange(ranges.object[i, ], bsSet)
	}
	state
}



onlyInOne <- function(distance.ranges, penn.ranges, CHR, xlim, method=c("d", "p"),
		      deletion=TRUE){
	if(deletion){
		states <- deletionStates()
	} else states <- duplicationStates()
	distance.ranges <- distance.ranges[distance.ranges$state %in% states, ]
	penn.ranges <- penn.ranges[penn.ranges$state %in% states, ]
	dist <- rangesInXlim(distance.ranges, CHR=CHR, xlim*1e6)
	penn <- rangesInXlim(penn.ranges, CHR=CHR, xlim*1e6)
	notInPenn <- dist$id[!dist$id %in% penn$id]
	notInDist <- penn$id[!penn$id %in% dist$id]
	if(method=="d") return(dist[dist$id %in% notInPenn, ])
	if(method=="p") return(penn[penn$id %in% notInDist, ])
	stop("method not 'p' or 'd'")
}


getTrioIndex <- function(bsSet){
	family=unique(sssampleNames(bsSet))
	family <- family[family!="CIDR"]
	family.wga <- sssampleNames(bsSet)[bsSet$is.wga]
	family <- family[!family %in% family.wga]
	f.index=grep("_03", sampleNames(bsSet))
	f.sns <- sssampleNames(bsSet)[f.index]
	##family <- family[match(f.sns, family)]
	##stopifnot(identical(f.sns, family))
	m.index=grep("_02", sampleNames(bsSet))
	m.sns <- sssampleNames(bsSet)[m.index]
	##family <- family[match(m.sns, family)]
	##stopifnot(identical(m.sns, family))
	o.index=grep("_01", sampleNames(bsSet))
	o.sns <- substr(sampleNames(bsSet)[o.index], 1, 5)
	##
	family.all <- intersect(intersect(o.sns, f.sns), m.sns)
	##stopifnot(identical(o.sns, family))
	f.index <- match(paste(family.all, "03", sep="_"), ssampleNames(bsSet))
	m.index <- match(paste(family.all, "02", sep="_"), ssampleNames(bsSet))
	o.index <- match(paste(family.all, "01", sep="_"), ssampleNames(bsSet))
	index <- cbind(f.index, m.index, o.index)
	colnames(index) <- c("F", "M", "O")
	return(index)
}

initTrioList <- function(bsSet, minDistanceSet){
	trioSet <- list(bsSet=bsSet,
			minDistanceSet=minDistanceSet,
			index=getTrioIndex(bsSet))
}

constructTrioSetFrom <- function(ranged.data,
				 trioList,
				 as.data.frame=TRUE,
				 FRAME,
				 ylim,
				 index.in.chromosome=TRUE,
				 verbose=TRUE,
				 unit="Mb",
				 ...){
	require(SNPchip)
	data(chromosomeAnnotation)
	open(logR(trioList[[1]]))
	open(baf(trioList[[1]]))
	## update the start/end indices with FRAME
	stopifnot(nrow(ranged.data)==1)
	##marker.index <- featuresInRange(trioList[["bsSet"]], ranged.data, )
	nms <- names(list(...))
	if("xlim" %in% nms){
		xlim <- list(...)[["xlim"]]
		if(unit=="Mb")
			xlim <- xlim*1e6
		rd2 <- addIndicesForRanges(ranged.data, trioList[["bsSet"]],
					   xlim=xlim,
					   index.in.chromosome=index.in.chromosome,
					   verbose=verbose)
	} else{
		rd2 <- addIndicesForRanges(ranged.data, trioList[["bsSet"]],
					   FRAME=FRAME,
					   index.in.chromosome=index.in.chromosome,
					   verbose=verbose)
	}
	if(index.in.chromosome){
		i <- rd2$start.index:rd2$end.index
		chr.index <- which(chromosome(trioList[["bsSet"]])==ranged.data$chrom)
		marker.index <- chr.index[i]
	}
	bsSet <- trioList$bsSet
	minDistanceSet <- trioList$minDistanceSet
	open(copyNumber(minDistanceSet))

	sampleNames(minDistanceSet) <- substr(sampleNames(minDistanceSet),1,5)
	j <- match(substr(ranged.data$id, 1,5), sampleNames(minDistanceSet))
	fmoindex <- match(paste(ss(ranged.data$id), c("03", "02", "01"), sep="_"), ssampleNames(bsSet))
	father.index <- fmoindex[1]##trioList$index[j, 1]
	mother.index <- fmoindex[2]##trioList$index[j, 2]
	offspring.index <- fmoindex[3]##trioList$index[j, 3]
	x <- list(
		   logR.F=as.matrix(logR(bsSet)[marker.index, father.index, drop=FALSE]),
		   logR.M=as.matrix(logR(bsSet)[marker.index, mother.index, drop=FALSE]),
		   logR.O=as.matrix(logR(bsSet)[marker.index, offspring.index, drop=FALSE]),
		   baf.F=as.matrix(baf(bsSet)[marker.index, father.index, drop=FALSE]),
		   baf.M=as.matrix(baf(bsSet)[marker.index, mother.index, drop=FALSE]),
		   baf.O=as.matrix(baf(bsSet)[marker.index, offspring.index, drop=FALSE]),
		   dist=as.matrix(copyNumber(minDistanceSet)[marker.index, j, drop=FALSE]))
	x <- lapply(x, function(x, sns) {colnames(x) <- sns; x}, sns=sampleNames(minDistanceSet)[j])
	open(minDistanceSet$MAD)
	obj <- new("MultiSet",
		   logR.F=x[["logR.F"]],
		   logR.M=x[["logR.M"]],
		   logR.O=x[["logR.O"]],
		   baf.F=x[["baf.F"]],
		   baf.M=x[["baf.M"]],
		   baf.O=x[["baf.O"]],
		   mindist=x[["dist"]],
		   phenoData=phenoData(minDistanceSet)[j, ],
		   featureData=featureData(bsSet)[marker.index, ])
	close(minDistanceSet$MAD)
	if(as.data.frame) obj <- as(obj, "data.frame")
	close(logR(trioList[[1]]))
	close(baf(trioList[[1]]))
	close(copyNumber(minDistanceSet))
	obj
}

findPeaks <- function(ranges.object){
	chrom <- unique(ranges.object$chrom)
	r.cnt <- vector("list", length(chrom))
	for(i in seq_along(chrom)){
		CHR <- chrom[i]
		r <- ranges.object[ranges.object$chrom==CHR,]
		r <- IRanges(start(r), end(r))
		dr <- disjointRanges(r)
		cnt <- countOverlaps(dr, r)
		r.cnt[[i]] <- RangedData(IRanges(start(dr), end(dr)),
					 freq=cnt,
					 chrom=CHR,
					 space=CHR)
	}
	do.call("c", r.cnt)
}

madVsCoverage <- function(lambda=0.1, MIN=1, MAX=4, coverage=3:100){
	p <- lambda*exp(-lambda*coverage) ## 0 - 0.04 (Pr (X=x)
	b <- 1/(MAX - MIN)
	a <- MIN * b
	numberMads <- ((p-min(p))/(max(p)-min(p)) + a)/b
	list(x=coverage, y=numberMads)
}

##calculateMinimumDistance <- function(trios, bsSet, minDistanceSet){


callTrioState <- function(ranges.object, subj.stat=c("F.median", "M.median", "O.median")){
	stopifnot(all(subj.stat %in% colnames(ranges.object)))
	triostate <- matrix(NA, nrow(ranges.object), 3)
	for(i in seq_along(subj.stat)){
		stat <- subj.stat[i]
		meds <- eval(substitute(ranges.object$NAME_ARG, list(NAME_ARG=stat)))
		homo.deletion <- meds <= -1
		hemi.deletion <- meds <= -0.1 & meds > -1
		normal <- meds > -0.1 & meds < 0.1
		amp <- meds >= 0.1
		tmp <- cbind(homo.deletion, hemi.deletion, normal, amp)
		which.state <- apply(tmp, 1, which)
		triostate[, i] <- as.integer(which.state)
	}
	triostate.i <- apply(triostate, 1, function(x) paste(x, collapse=""))
	as.integer(triostate.i)

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

mset2df <- function(mset, ylim){
	df <- as(mset, "data.frame")
	index <- which(df$logR < ylim[1])
	df$logR[index] <- jitter(rep(ylim[1], length(index)), amount=0.1)
	df
}

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

distanceRangesForRange <- function(range.object, distance.ranges){
	stopifnot(nrow(range.object)==1)
	CHR <- range.object$chrom[[1]]
	id <- substr(range.object$id, 1, 5)
	id2 <- substr(distance.ranges$id, 1, 5)
	index <- id2 %in% id & distance.ranges$chrom == CHR
	distance.ranges[index, ]
}

dropHemizygousRangesWithHets <- function(range.object, pHet.thr=0.01){
	stopifnot("triostate" %in% colnames(range.object))
	hemizygous.index <- which(substr(range.object$triostate, 3, 3)==2)
	het.index <- which(range.object$pHet > pHet.thr)
	drop.index <- intersect(hemizygous.index, het.index)
	range.object <- range.object[-drop.index, ]
	range.object
}

dropDuplicatedRangesWithHets <- function(range.object, pHet.thr=0.01, bsSet){
	stopifnot("triostate" %in% colnames(range.object))
	duplicated.index <- which(substr(range.object$triostate, 3, 3) > 3)
	range.object <- range.object[duplicated.index, ]
	mn <- rep(NA, length(duplicated.index))
	for(i in seq_along(duplicated.index)){
		j <- duplicated.index[i]
		marker.index <- (range.object$start.index[j]):(range.object$end.index[j])
		sample.index <- match(range.object$id[j], sampleNames(bsSet))
		b <- as.matrix(baf(bsSet)[marker.index, sample.index])
		mn[i] <- mean(b > 0.45 & b < 0.55)
	}
	het.index <- which(range.object$pHet > pHet.thr)
	drop.index <- intersect(hemizygous.index, het.index)
	range.object <- range.object[-drop.index, ]
	range.object
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
##		if(any(st >= en)) {
##			index <- which(st < en & end(range.object) > en)
##			range.object <- range.object[index, ]
##			en <- en[index]
##			st <- st[index]
##		}
##		stopifnot(st < en)
##		end(range.object) <- en
	}
	if(insertWhere=="centromereEnd"){
		##st <- chromosomeAnnotation[chrom, j]
		##en <- end(range.object)
		##if(any(st >= en)) {
		index <- which(start(range.object) <= chromosomeAnnotation[chrom, 2] & end(range.object) >= chromosomeAnnotation[chrom, 2])
		stopifnot(length(index) == 1)
		##which(chromosomeAnnotation[chrom, j] < end(range.object) & chromosomeAnnotation[chrom, j]
		##index <- which(st < en & start(range.object) >= st)
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

##	rbind(range.object, ranges.updated)
##}
departureFromZero <- function(minDistanceSet,
			      minDistanceRanges,
			      bsSet,
			      MAX=4,
			      MIN=1.5,
			      MIN.THR=0.25,
			      MIN.COVERAGE=5,
			      lambda=0.1,
			      pHet.thr=0.01){
##	xy <- madVsCoverage(lambda, MIN, MAX)
##	xyplot(y~x, xy, type="l", ylab="# of MADs", xlab="coverage", xlim=c(0,60),
##	       scales=list(num.tick=10))
##	trellis.focus("panel", 1, 1)
##	panel.abline(v=c(3, 10, 20, 30), lty=3, col="grey")
##	deletion.ranges <- collectAllRangesOfSize(SIZE=MIN.COVERAGE,
##						  bsSet=bsSet,
##						  minDistanceSet=minDistanceSet,
##						  minDistanceRanges=minDistanceRanges,
##						  MIN=MIN,
##						  MAX=MAX,
##						  lambda=lambda,
##						  denovo.deletion=TRUE,
##						  MIN.THR=MIN.THR)
##	deletion.ranges$type <- "denovo_deletion"
##	amp.ranges <- collectAllRangesOfSize(SIZE=MIN.COVERAGE,
##					     bsSet=bsSet,
##					     minDistanceSet=minDistanceSet,
##					     minDistanceRanges=minDistanceRanges,
##					     MIN=MIN,
##					     MAX=MAX,
##					     lambda=lambda,
##					     MIN.THR=MIN.THR,
##					     denovo.deletion=FALSE)
##	amp.ranges$type <- "denovo_duplication"
	## drop all ranges that are not 1.5 mads from zero
	open(minDistanceSet$MAD)
	i <- match(minDistanceRanges$id, sampleNames(minDistanceSet))
	mads <- minDistanceSet$MAD[i]
	cutoff <- 2*mads
	is.altered <- abs(minDistanceRanges$seg.mean) > cutoff
	altered.ranges <- minDistanceRanges[is.altered, ]
	altered.ranges <- altered.ranges[altered.ranges$num.mark >= MIN.COVERAGE, ]
	altered.ranges <- replaceCentromereRanges(altered.ranges)
	##post-hoc calls
	message("computing post-hoc calls for the trio state...")
	triostate <- callTrioState(altered.ranges)
	altered.ranges$triostate <- triostate
	##altered.ranges <- altered.ranges[triostate %in% c(deletionStates(), duplicationStates()), ]
	## filter hets in hemizygous deletions
	message("dropping chemizygous offspring calls with prop(AB) exceeding ", pHet.thr)
	altered.ranges <- dropHemizygousRangesWithHets(altered.ranges, pHet.thr)
	##altered.ranges <- dropDuplicatedRangesWithHets(altered.ranges, pHet.thr)
	return(altered.ranges)
}

initializeffForEmission <- function(outdir, nr, nc){
	logEP.cn1 <- initializeBigMatrix("emission.cn1", nr=nr, nc, vmode="double")
	logEP.cn2 <- initializeBigMatrix("emission.cn2", nr=nr, nc, vmode="double")
	logEP.cn3 <- initializeBigMatrix("emission.cn3", nr=nr, nc, vmode="double")
	logEP.cn4 <- initializeBigMatrix("emission.cn4", nr=nr, nc, vmode="double")
	logEP.cn <- ffdf(logEP.cn1, logEP.cn2, logEP.cn3, logEP.cn4)

	logEP.b1 <- initializeBigMatrix("emission.b1", nr=nr, nc, vmode="double")
	logEP.b2 <- initializeBigMatrix("emission.b2", nr=nr, nc, vmode="double")
	logEP.b3 <- initializeBigMatrix("emission.b3", nr=nr, nc, vmode="double")
	logEP.b4 <- initializeBigMatrix("emission.b4", nr=nr, nc, vmode="double")
	logEP.bf <- ffdf(logEP.b1, logEP.b2, logEP.b3, logEP.b4)
	list(copynumber=logEP.cn, baf=logEP.bf)
}

musInRange <- function(query, cbs.segs, id, chr){
	index <- which(cbs.segs$id == id)
	cbs.sub <- cbs.segs[index, ]
	subj <- IRanges(start(cbs.sub), end(cbs.sub))
	mm <- matchMatrix(findOverlaps(query, subj))[,2]
	cbs.sub[mm, ]
}

frameRange <- function(range.object){

}

intervalOfAlterations <- function(altered.ranges, bsSet){
	I <- nrow(bsSet)
	start.index <- altered.ranges$start.index
	end.index <- altered.ranges$end.index
	CHR <- chromosome(bsSet)[start.index]
	previousI <- start.index - 1
	nextI <- end.index + 1
	chromPreviousI <- chromosome(bsSet)[previousI]
	chromNextI <- chromosome(bsSet)[nextI]
	previousI <- ifelse(chromPreviousI != CHR, start.index, previousI)
	nextI <- ifelse(chromNextI != CHR, end.index, nextI)
	indexInterval <- cbind(previousI, start.index, end.index, nextI, CHR)
	colnames(indexInterval) <- c("previous.index", "start.index", "end.index", "next.index", "chrom")
	indexInterval
##	posRange <- cbind(position(bsSet)[indexInterval[, 1]], position(bsSet)[indexInterval[,2]])
##	posRange
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

calculateChangeSd <- function(coverage=1:500, lambda, a=0.2, b=0.025)
	a + lambda*exp(-lambda*coverage)/b


prune <- function(genomdat,
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
	trimmed.SD <- unique(range.object$mad)
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
			##
			##thrSD <- y*trimmed.SD
			##thrSD[thrSD < MIN.CHANGE] <- MIN.CHANGE
			segments0 <- cbind(c(1,1+cpt.loc[-k]),cpt.loc)
			## now:
			## 1   143
			## 144 152
			## 153 165
			## median copy number for each segment
			segmed <- apply(segments0, 1, function(i,x) {median(x[i[1]:i[2]], na.rm=T)}, genomdat)
			## coverage <- apply(segments0, 1, function(i) length(i[1]:i[2]))
			## coverage is the same as cpt.loc...unnecessary
			##
			## absolute copy number difference of adjacent segments
 			##adsegmed <- abs(diff(segmed))
			adsegmed <- abs(diff(segmed))  ## put abs in inner parentheses to keep things far from zero
			##
			##
			## number of standard deviations of observed shift
			empiricalNumberSd <- adsegmed/trimmed.SD
			##
			## add to this the difference of the last segment and zero.
			##stopifnot(identical(length(adsegmed), length(thrSD)))
			## drop order: coverage then distance
			## i <- order(coverage, adsegmed)[1]
##			if(any(adsegmed < thrSD | coverage < MIN.COVERAGE)){
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
		   id=unique(range.object$id),
		   chrom=unique(range.object$chrom),
		   num.mark=lseg,
		   seg.mean=segmeans,
		   start.index=segments0[,1],
		   end.index=segments0[,2],
		   space=id)
}

pruneRanges <- function(CHR,
			bsSet,
			minDistanceSet,
			ids,
			lambda=0.05,
			MIN.CHANGE=0.1,
			SCALE.EXP=0.02,
			MIN.COVERAGE=3, verbose=TRUE){
	stopifnot(length(CHR)==1)
	open(copyNumber(minDistanceSet))
	open(baf(bsSet))
	open(logR(bsSet))
	open(bsSet$MAD)
	open(minDistanceSet$MAD)
	range.object <- loadRanges(beadstudiodir(), pattern=paste("md.segs_chr", CHR, "_", sep=""), CHR=CHR, name="md.segs")
	if(!missing(ids)) range.object <- range.object[ss(range.object$id) %in% ids, ]
	fD <- featureData(bsSet)[chromosome(bsSet)==CHR, ]
	range.object <- addIndicesFromFeatureData(range.object, fD)
	mads <- minDistanceSet$MAD
	mads <- mads[match(range.object$id, sampleNames(minDistanceSet))]
	range.object$mad <- mads
	## pruning is sample-specific (different noise)
	samples.in.range <- unique(range.object$id)
	rdList <- vector("list", length(samples.in.range))
	for(j in seq_along(samples.in.range)){
		if(verbose) if(j %% 100==0) message('   sample (', j, '/', length(samples.in.range), ')')
		rd <- range.object[range.object$id==samples.in.range[j], ]
		k <- match(samples.in.range[j], sampleNames(minDistanceSet))
		genomdat <- as.numeric(copyNumber(minDistanceSet)[chromosome(minDistanceSet) == CHR, k])
		##trace(prune, browser)
		range.pruned <- prune(genomdat,
				      rd,
				      physical.pos=fD$position,
				      lambda=lambda,
				      MIN.CHANGE=MIN.CHANGE,
				      SCALE.EXP=SCALE.EXP,
				      MIN.COVERAGE=MIN.COVERAGE)
		rdList[[j]] <- range.pruned
	}
	rdList <- rdList[!sapply(rdList, is.null)]
	range.pruned <- do.call("c", rdList)
	range.pruned <- RangedData(IRanges(start(range.pruned), end(range.pruned)),
				    chrom=range.pruned$chrom,
				    id=range.pruned$id,
				    num.mark=range.pruned$num.mark,
				    seg.mean=range.pruned$seg.mean,
				    start.index=range.pruned$start.index,
				    end.index=range.pruned$end.index)
	close(copyNumber(minDistanceSet))
	close(baf(bsSet))
	close(logR(bsSet))
	close(bsSet$MAD)
	close(minDistanceSet$MAD)
	return(range.pruned)
}

##prune <- function(range.object, minDistanceSet, ...){
##	drop.col <- match("pedId", colnames(range.object))
##	if(!is.na(drop.col))
##		range.object <- range.object[, -drop.col]
##	CHR <- unique(range.object$chrom)
##	stopifnot(length(CHR) == 1)
##	fD <- featureData(bsSet)[chromosome(minDistanceSet)==CHR, ]
##	##trace(addIndicesFromFeatureData, browser)
##	message("Adding indices from chromosome ", CHR, " featureData...")
##	range.object <- addIndicesFromFeatureData(range.object, fD)
##	sns <- sampleNames(minDistanceSet)
##	rdList <- vector("list", ncol(minDistanceSet))
##	for(j in 1:ncol(minDistanceSet)){
##		if(j %% 100 == 0) cat(".")
##		id <- sns[j]
##		index.ranges <- which(range.object$id == id)
##		if(length(index.ranges) > 1){
##			##stop()
##			genomdat <- as.numeric(copyNumber(minDistanceSet)[chromosome(minDistanceSet) == CHR, j])
##			##trace(myprune, browser)
##			rdList[[j]] <- myprune(genomdat, range.object[index.ranges, ], physical.pos=fD$position, trimmed.SD=mads[j],...)
##		} else rdList[[j]] <- range.object
##	}
##	range.pruned <- do.call("c", rdList)
##	range.pruned
##}
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
phi <- function(x, mu, sigma) dnorm(x, mu, sigma)
## cdf of standard normal
Phi <- function(x, mu, sigma) pnorm(x, mu, sigma)
## pdf of truncated normal on support [0, 1]
tnorm <- function(x, mu, sigma) phi(x, mu, sigma)/(Phi(1, mu, sigma)-Phi(0, mu, sigma))


likelihood <- function(CHR, bsSet, index.list, family.list, chr.index, states=0:4){
	lname <- file.path(ldPath(), paste("L_", CHR, ".rda", sep=""))
	if(!file.exists(lname)){
		message("Initializing new ff objects in ", ldPath())
		LF <- createFF(paste("likF_chr", CHR, "_", sep=""),
			       dim=c(length(unlist(index.list)),
			       length(unlist(family.list)),
			       length(states)), vmode="double", initdata=NA)
		LM <- createFF(paste("likM_chr", CHR, "_", sep=""),
			       dim=c(length(unlist(index.list)),
			       length(unlist(family.list)),
			       length(states)), vmode="double", initdata=NA)
		LO <- createFF(paste("likO_chr", CHR, "_", sep=""),
			       dim=c(length(unlist(index.list)),
			       length(unlist(family.list)),
			       length(states)), vmode="double", initdata=NA)
		colnames(LF) <- colnames(LM) <- colnames(LO) <- unlist(family.list)
		L <- list(F=LF, M=LM, O=LO)
		save("L", file=lname)
	} else{
		message("Loading ", lname)
		load(lname)
		LF <- L[["F"]]
		LM <- L[["M"]]
		LO <- L[["O"]]
	}
	for(j in seq_along(family.list)){
		##J <- which(sssampleNames(bsSet) %in% family.list[[j]])
		fams <- family.list[[j]]
		##i.fam <- match(paste(fams, c("03", "02", "01"), sep="_"), ssampleNames(bsSet))
		i.f <- match(paste(fams, "03", sep="_"), ssampleNames(bsSet))
		ocLapply(seq_along(index.list),
			 calculateLikelihood,
			 bsSet=bsSet,
			 L=LF,
			 index.list=index.list,
			 j=i.f,
			 chr.index=chr.index,
			 states=states)
		i.m <- match(paste(fams, "_02", sep=""), ssampleNames(bsSet))
		ocLapply(seq_along(index.list),
			 calculateLikelihood,
			 bsSet=bsSet,
			 L=LM,
			 index.list=index.list,
			 j=i.m,
			 chr.index=chr.index,
			 states=states)
		i.o <- match(paste(fams, "_01", sep=""), ssampleNames(bsSet))
		ocLapply(seq_along(index.list),
			 calculateLikelihood,
			 bsSet=bsSet,
			 L=LO,
			 index.list=index.list,
			 j=i.o,
			 chr.index=chr.index,
			 states=states)
	}
	list(F=LF,
	     M=LM,
	     O=LO)
}

constructSet <- function(bsSet, CHR, id, states){
	open(baf(bsSet))
	open(logR(bsSet))
	open(bsSet$MAD)
	i <- which(chromosome(bsSet) == CHR)
	j <- match(id, ssampleNames(bsSet))
	S <- length(states)
	loglik <- array(NA, dim=c(2, length(i), length(j), S))
	dimnames(loglik) <- list(c("logR", "baf"),
			      featureNames(bsSet)[i],
			      sampleNames(bsSet)[j],
			      states)
	object <- new("LikSet",
		      logR=as.matrix(logR(bsSet)[i,j]),
		      BAF=as.matrix(baf(bsSet)[i,j]),
		      phenoData=phenoData(bsSet)[j, ],
		      featureData=featureData(bsSet)[i, ],
		      experimentData=experimentData(bsSet),
		      annotation=annotation(bsSet),
		      protocolData=protocolData(bsSet)[j, ],
		      loglik=loglik)
	fData(object)$range.index <- NA
	close(baf(bsSet))
	close(logR(bsSet))
	close(bsSet$MAD)
	return(object)
}

computeLoglik <- function(id,
			  bsSet, #L,
			  CHR,
			  mu.logr=c(-2, -0.5, 0, 0.3, 0.75),
			  states=0:4,
			  baf.sds=c(0.02, 0.03, 0.02),
			  THR=-50,
			  prGtCorrect=0.999){ ##prob genotype is correct
	p1 <- prGtCorrect; rm(prGtCorrect)
	##
	## one obvious thing that p1 could depend on is the
	## minor allele frequency.  If rare, p1 is smaller
	##
	##i <- index.list[[stratum]]
	##k <- chr.index[[stratum]]
	stopifnot(all(!is.na(match(id, ssampleNames(bsSet)))))
	object <- constructSet(bsSet, CHR, id, states=states)
##	lR <- as.matrix(logR(bsSet)[i, j])
##	bf <- as.matrix(baf(bsSet)[i, j])
	##lik.logr <- array(NA, dim=c(nrow(object), ncol(object), length(states)))
	##lik.baf <- array(NA, dim=dim(lik.logr))
##	lik <- array(NA, dim=c(2, nrow(object), ncol(object), length(states)))
##	dimnames(lik) <- list(c("logR", "baf"),
##			      featureNames(object),
##			      sampleNames(object),
##			      states)
	open(bsSet$MAD)
	med.all <- median(bsSet$MAD[!bsSet$is.wga], na.rm=TRUE)
	j <- match(id, ssampleNames(bsSet))
	sds.sample <- bsSet$MAD[j]
##	scale.sd <- sds.sample
##	scale.sd[sds.sample < med.all] <- 1
##	scale.sd[sds.sample >= med.all] <- sds.sample[sds.sample >= med.all]/med.all
	##scale.sd <- ifelse(sds.sample < median(sds.sample), 1, sds.sample/median(sds.sample))
##	scale.sd <- matrix(scale.sd, nrow(object), length(j), byrow=TRUE)
	##sds.marker <- crlmm:::rowMAD(lR)
	##marker.index <- which(chromosome(bsSet) == CHR)
	sds.sample <- matrix(sds.sample, nrow(object), length(j), byrow=TRUE)
	sds.marker <- fData(object)$MAD
	sds.marker <- matrix(sds.marker, nrow(object), length(j), byrow=FALSE)
	## less ad hoc: shrink the marker-specific estimate to the sample-level estimate
	sds <- (sds.marker+sds.sample)/2
	lR <- logR(object)
	##lik <- loglik(object)
	for(i in seq_along(states)) loglik(object)["logR", , , i] <- dnorm(lR, mu.logr[i], sds)
##	loglik(object)["logR", 12572:12590, 1, c(2,3)]
##	s1 <- apply(loglik(object)["logR", 12572:12590, 1, c(2,3)], 2, sum)
##	s2 <- apply(loglik(object)["logR", 12572:12590, 2, c(2,3)], 2, sum)
##	s3 <- apply(loglik(object)["baf", 12572:12590, 1, c(2,3)], 2, sum)
##	s4 <- apply(loglik(object)["baf", 12572:12590, 2, c(2,3)], 2, sum)
##	lik["logR", , 1] <- dnorm(lR, -2, sds)
##	lik["logR", , 2] <- dnorm(lR, -0.5, sds)
##	lik["logR", , 3] <- dnorm(lR, 0, sds)
##	lik["logR", , 4] <- dnorm(lR, 0.3, sds)
##	lik["logR", , 5] <- dnorm(lR, 0.75, sds)
	sd0 <- baf.sds[1]
	sd.5 <- baf.sds[2]
	sd1 <- baf.sds[3]
	## we can also multiply the density given by tnorm by the
	## binomial probability of the genotype as per Wang.
	##  -- could use estimates of the allele frequency for the p parameter
	bf <- baf(object)
	## model emission as a mixture of normals (genotype is correct) and a uniform (error)
	perr <- (1-p1)*runif(bf)
	## Wang et al. use mixture probabilities from a binomial
	##
	##   pi ~ binomial(C(z), allele freq)
	##   and integrates out the number of copies of the B allele.
	##
	##   Below, I,ve just used a mixture model.  I have not integrated out the	#      copy number of the B allele, nor do I make use of MAF estimates.
	loglik(object)["baf", , , 1] <-  1
	loglik(object)["baf", , , 2] <- p1*(1/2*tnorm(bf, 0, sd0) + 1/2*tnorm(bf, 1, sd1)) + perr
	loglik(object)["baf", , , 3] <- p1*(1/3*tnorm(bf, 0, sd0) + 1/3*tnorm(bf, 0.5, sd.5) + 1/3*tnorm(bf, 1, sd1))+perr
	loglik(object)["baf", , , 4] <- p1*(1/4*tnorm(bf, 0, sd0) + 1/4*tnorm(bf, 1/3, sd.5) + 1/4*tnorm(bf, 2/3, sd.5) + 1/4*tnorm(bf, 1, sd1)) + perr
	loglik(object)["baf", , , 5] <- p1*(1/5*tnorm(bf, 0, sd0) + 1/5*tnorm(bf, 1/4, sd.5) + 1/5*tnorm(bf, 0.5, sd.5) + 1/5*tnorm(bf, 3/4, sd.5) + 1/5*tnorm(bf, 1, sd1)) + perr
	loglik(object) <- log(loglik(object))
	loglik(object)[loglik(object) < THR] <- THR
	##ll[ll < THR] <- THR
	##loglik(object) <- ll
	## assume an epsilon probability that the observation was an error
	##     (truncating the lik. should have a similar effect)
##	loglik[loglik < THR] <- THR
	##close(L)
	close(bsSet$MAD)
	return(object)
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

##transitionProbability <- function(states, epsilon=1-0.999){
##	off.diag <- epsilon/(nrow(states) -1)
##	tpm <- matrix(off.diag, nrow(states), nrow(states))
##	diag(tpm) <- 1-epsilon
##	tpm
##}

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

joint1 <- function(LLT, ##object,
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
	##browser()
	state <- trio.states[state.index, ]
	##fmo <- list()
	##fmo <- matrix(NA, nrow(object), 3)
	##tmp <- as.matrix(do.call("cbind", fmo))
	##for(i in 1:3) fmo[, i] <- loglik(object)["logR", , i, state[i]] + loglik(object)["baf", , i, state[i]]
	##fmo <- fmo[rowSums(is.na(fmo)) == 0, ]
	##LLT is 3 x 5
	fmo <- c(LLT[1, state[1]], LLT[2, state[2]], LLT[3, state[3]])
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
		##fmo <- apply(fmo, 2, sum, na.rm=TRUE)
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
		##fmo <- apply(fmo, 2, sum, na.rm=TRUE)
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
	## RS: comment 4/29/2011
	## add the probability of transitioning back to the normal state
	##for(j in 1:3) fmo[j] <- fmo[j] + log(tau[state[j], 3])
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

joint4 <- function(bsSet,
		   ranges,
		   states,
		   baf.sds,
		   THR=-50,
		   mu.logr=c(-2,-0.5, 0, 0.3, 0.75),
		   log.pi,
		   tau,
		   normal.index,
		   a=0.0009,
		   verbose=TRUE,
		   prGtCorrect=0.999){
	stopifnot(states == 0:4)
	##ids <- unique(ss(ranges$id))
	##sample.id <- ids[family.index]
	##sample.index <- match(sample.id, colnames(F))
	family.id <- unique(ss(ranges$id))
	fmonames <- paste(ss(family.id), c("03", "02", "01"), sep="_")
	object <- computeLoglik(id=fmonames,
				bsSet=bsSet,
				CHR=ranges$chrom[[1]],
				mu.logr=mu.logr,
				states=states,
				baf.sds=baf.sds,
				THR=THR,
				prGtCorrect=prGtCorrect)
	##ll <- loglik(object)
	start.stop <- cbind(ranges$start.index, ranges$end.index)
	l <- apply(start.stop, 1, function(x) length(x[1]:x[2]))
	fData(object)$range.index <- rep(seq(length=nrow(ranges)), l)
	##rtmp <- loglik(object)["logR", range.index(object)==2, ,  ]
	##btmp <- loglik(object)["baf", range.index(object)==2, ,  ]
	trio.states <- trioStates(states)
	tmp <- matrix(NA, nrow(trio.states), 2)
	colnames(tmp) <- c("DN=0", "DN=1")
	## take into account 'the prior'
	##  -- the probability of transitioning to and from an altered state for each range
	##  -- the initial state probability if range is denovo
	##  -- the initial state probability of range is not denovo
	state.prev <- NULL
	denovo.prev <- NULL
	table1 <- readTable1(a=a)
	table3 <- readTable3(a=a)
	weightR <- 1/3
	for(i in seq(length=nrow(ranges))){
		obj <- object[range.index(object) == i, ]
		LLR <- loglik(obj)["logR", , ,  ]
		LLB <- loglik(obj)["baf", , , ]
		LL <- weightR * LLR + (1-weightR)*LLB
		##
##		## recenter and scale each
##		LLR.c <- LLR
##		LLB.c <- LLB
##		for(j in 1:3){
##			## centering and scaling the rows gives each probe equal weight
##			## gives the approx. BAF equal importance
##			x <- t(LLR[, j, ])
##			tx <- t(scale(x))
##			LLR.c[, j, ] <- tx
##			x <- t(LLB[, j, ])
##			tx <- t(scale(x))
##			LLB.c[, j, ] <- tx
##		}
##		LL.c <- LLR.c+LLB.c
		##
		## Does it make sense to give the logR and BAFs emission probs equal weight?
		## - The logR emission probs are much more variable and can dominate
		##
		##
##		LL <-  LLR + LLB
##		LL[is.na(LL)] <- log(1e-10)
		## nr x 3, S
		## sum over rows -> 3 x S
		LLT <- matrix(NA, 3, 5)
		for(j in 1:3) LLT[j, ] <- apply(LL[, j, ], 2, sum, na.rm=TRUE)
		for(j in 1:nrow(trio.states)){
			for(DN in c(FALSE, TRUE)){
				tmp[j, DN+1] <- joint1(##object=obj,
						       LLT=LLT,
						       trio.states=trio.states,
						      tau=tau,
						      log.pi=log.pi,
						      normal.index=normal.index,
						      segment.index=i,
						      state.index=j,
						      table1=table1,
						      table3=table3,
						      is.denovo=DN,
						      state.prev=state.prev,
						      denovo.prev=denovo.prev)
			}
		}
		## RS 4/29/2011
		##integrate out the denovo indicator
		##one.finite <- which(rowSums(is.finite(tmp))==1)
		argmax1 <- which.max(tmp[,1])
		argmax2 <- which.max(tmp[,2])
		if(argmax1 != argmax2){
			lik1 <- tmp[argmax1, 1]
			lik2 <- tmp[argmax2, 2]
			if(lik1 >= lik2){
				argmax <- argmax1
				is.denovo <- FALSE
				bf <- tmp[argmax1, 1]
			} else{
				is.denovo <- TRUE
				argmax <- argmax2
				bf <- tmp[argmax2, 2]
			}
		} else{
			argmax <- argmax1
			is.denovo <- FALSE
			bf <- tmp[argmax1, 1]
		}
##		##tmp[!is.finite(tmp)] <- min(tmp[is.finite(tmp)], na.rm=TRUE)
##		## in most cases, a non-mendelian event should always be finite
##		##  -- only the mendelian events will be -Inf
##		##     -Inf + finite = -Inf
##		##     (this will automatically exclude a state when a mendelian mechanism is improbable)
##		##     -> -Inf should only have an effect when the true mechanism is non-mendelian
##		##     -> If in truth the mechanism is non-mendelian, the right answer would be the max of column 2
##		##     -> If using rowSums,  we're choosing the state that could arise by a mendelian mechanism
##		rsums <- rowSums(tmp, na.rm=TRUE)
##		if(all(is.infinite(rsums))){
##			min.val <- min(tmp[is.finite(tmp)], na.rm=TRUE)
##			tmp[is.infinite(tmp)] <- min.val
##			rsums <- rowSums(tmp, na.rm=TRUE)
##			##stop("all -Inf in joint4")
##		}
##		##argmax1 <- which.max(tmp[, 1])
##		##argmax2 <- which.max(tmp[, 2])
##		##argmax <- which.max(tmp)
##		##if(is.denovo){
##		##bf <- tmp[argmax2, 2]-tmp[normal.index, 1]
##		##argmax <- argmax2
##		##}  else {
##		argmax <- which.max(rsums)
##		is.denovo <- ifelse(tmp[argmax, 1] < tmp[argmax, 2], TRUE, FALSE)
##		bf <- rsums[argmax] - rsums[normal.index]
		##bf <- tmp[argmax1, 1]-tmp[normal.index, 1]
		##argmax <- argmax1
		ranges$bayes.factor[i] <- bf
		ranges$DN[i] <- is.denovo
		ranges$argmax[i] <- argmax
		denovo.prev <- is.denovo
		state.prev <- trio.states[argmax, ]
	}
	ranges
}

computeBayesFactor <- function(range.object,
			       bsSet,
			       states=0:4,
			       baf.sds=c(0.02, 0.03, 0.02),
			       THR=-50,
			       mu.logr=c(-2, -0.5, 0, 0.3, 0.75),
			       log.pi,
			       tau,
			       normal.index=61,
			       a=0.0009,
			       prGtCorrect=0.999,
			       verbose=TRUE){
	stopifnot(!missing(tau))
	stopifnot(!missing(log.pi))
	stopifnot(!missing(bsSet))
	range.object$bayes.factor <- NA
	range.object$argmax <- NA
	range.object$DN <- NA
	ids <- unique(range.object$id)
	for(i in seq_along(ids)){
		this.id <- ids[i]
		if(verbose){
			if(i %% 100 == 0)
				message("   sample ", this.id, " (", i, "/", length(ids), ")")
		}
		j <- which(range.object$id == this.id)
		rd <- joint4(bsSet=bsSet,
			     ranges=range.object[j, ],
			     states=states,
			     baf.sds=baf.sds,
			     THR=THR,
			     mu.logr=mu.logr,
			     log.pi=log.pi,
			     tau=tau,
			     normal.index=normal.index,
			     a=a,
			     prGtCorrect=prGtCorrect,
			     verbose=verbose)##, F=F, M=M, O=O)
		range.object$bayes.factor[j] <- rd$bayes.factor
		range.object$argmax[j] <- rd$argmax
		range.object$DN[j] <- rd$DN
	}
	range.object
}

pruneThenBayesFactor <- function(CHR,
				 bsSet,
				 minDistanceSet,
				 ids,
				 lambda=0.05,
				 MIN.CHANGE=0.1,
				 SCALE.EXP=0.02,
				 MIN.COVERAGE=10,
				 states=0:4,
				 baf.sds=c(0.02, 0.03, 0.02),
				 THR=-50,
				 mu.logr=c(-2, -0.5, 0, 0.3, 0.75),
				 log.pi,
				 tau,
				 normal.index=61,
				 a=0.0009,
				 prGtCorrect=0.999,
				 verbose=TRUE){
	stopifnot(!missing(CHR))
	##minDistanceRanges <- minDistanceRanges[minDistanceRanges$chrom == CHR, ]
	if(verbose) message("Pruning chromosome ", CHR)
	stopifnot(!missing(tau))
	range.pruned <- pruneRanges(CHR=CHR,
				    bsSet=bsSet,
				    minDistanceSet=minDistanceSet,
				    ids=ids,
				    lambda=lambda,
				    MIN.CHANGE=MIN.CHANGE,
				    SCALE.EXP=SCALE.EXP,
				    MIN.COVERAGE=MIN.COVERAGE,
				    verbose=verbose)
	if(verbose) message("Computing bayes factors for chromosome ", CHR)
	ranges.bf <- computeBayesFactor(range.object=range.pruned,
					bsSet=bsSet,
					states=states,
					baf.sds=baf.sds,
					THR=THR,
					mu.logr=mu.logr,
					log.pi=log.pi,
					tau=tau,
					normal.index=normal.index,
					a=a,
					prGtCorrect=prGtCorrect,
					verbose=verbose)
	## consecutive ranges that have the same state can be
	## collapsed into 1 range with the posterior odds added (or
	## multiplied)
	##
	##  -- a less ad-hoc approach is to collapse and recompute the
	##     bayes factor, but this is also more time-consuming.
	ranges.bf <- pruneByFactor(ranges.bf, f=ranges.bf$argmax)
	if(verbose) message("Collapsing adjacent ranges with same trio state")
	if(is(ranges.bf, "list")) {
		## each element should be RangedData
		tmp <- tryCatch(res <- RangedDataList(ranges.bf), error=function(e) NULL)
		##tmp <- tryCatch(res <- do.call("rbind", ranges.bf), error=function(e) NULL)
		if(!is.null(tmp)) {
			res <- stack(res)
		} else	res <- ranges.bf
	} else res <- ranges.bf
	return(res)
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
		sds <- crlmm:::rowMAD(lr, na.rm=TRUE)
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

p1 <- function(df, palette, xlim, show.coverage=TRUE, blackBorder=TRUE, sampleLabels.cex=0.5, labelAllSamples=TRUE,...){
	stopifnot(length(unique(df$chr))==1)
	mykey <- simpleKey(c("homo-del", "hemi-del", "amp")[palette %in% df$col], points=FALSE,
		   rectangles=TRUE, col=palette[palette %in% df$col], space="top")
	mykey$rectangles[["border"]] <- mykey$rectangles[["col"]] <- palette[palette %in% df$col]
	##df$method=factor(df$method, order=TRUE)
	if(blackBorder) border <- rep("black", nrow(df)) else border <- df$col
	if(labelAllSamples) {
		labels <- df$id
		ticks.at <- df$y
	} else {
		labels <- FALSE
		ticks.at <- pretty(df$y)
	}
	fig <- xyplot(y~midpoint|method, data=df,
			    panel=function(x, y, x0, x1, chr.size,
			    col, border, coverage, chr, show.coverage=TRUE, max.y,
			    ..., subscripts){
				    panel.grid(h=-1, v=10)
				    ##yy <- factor(y, order=TRUE)
				    panel.xyplot(x, y, ..., subscripts)
				    ##yyy <- as.integer(yy)
				    h <- 0.75
				    lrect(xleft=x0[subscripts],
					  xright=x1[subscripts],
					  ybottom=y-h/2,
					  ytop=y+h/2,
					  border=border[subscripts],
					  col=col[subscripts], ...)
				    if(show.coverage)
					    ltext(x, y,labels=coverage[subscripts], cex=0.6)
				    ## plot centromere
				    chr <- unique(as.integer(as.character(df$chr)))
				    coords <- chromosomeAnnotation[chr, 1:2]/1e6
				    lrect(xleft=coords[1],
					  xright=coords[2],
					  ybottom=0,
					  ytop=max.y+h/2,
					  col="grey",
					  border="grey")
			    },
		      x0=df$x0,
		      x1=df$x1,
		      col=df$col,
		      border=border,
		      alpha=1,
		      chr.size=df$chr.size,
		      scales=list(y=list(labels=labels, at=ticks.at, cex=sampleLabels.cex)),
		      coverage=df$coverage,
		      xlab="Mb",
		      ylab="offspring index",
		      show.coverage=show.coverage,
		      key=mykey,
		      par.strip.text=list(lines=0.7, cex=0.6),
		      prepanel=prepanel.fxn,
		      xlim=xlim,
		      max.y=max(df$y), ...)
##		      axis=function(side, text.cex){
##			      panel.axis(side, text.cex=text.cex)}, ...)
	return(fig)
}

logrpanelfunction2 <- function(x, y,
			       col,
			       fill,
			       ylimit,
			       highlight,
			       add.segments=TRUE,
			       range.object,
			       cbs.segs,
			       dranges,
			       alpha=0.5,
			       highlight.col=makeTransparent("blue",alpha=alpha),
			       ..., subscripts){
	stopifnot(nrow(range.object) == 1)
	if(missing(cbs.segs)){
		CHR <- range.object$chrom
		cbs.segs <- loadRangesCbs(beadstudiodir(), pattern=paste("cbs_chr", CHR, ".rda", sep=""), CHR=CHR, name="cbs.segs")
	}
 	stopifnot(nrow(cbs.segs)>0)
	if(missing(dranges)){
		dranges <- loadRanges(beadstudiodir(), pattern=paste("md.segs_chr", CHR, "_", sep=""), CHR=CHR, name="md.segs")
		dranges <- dranges[dranges$id==range.object$id[[1]], ]
	}
	stopifnot(nrow(dranges)>0)
	##stopifnot(nrow(cbs.segs) > 0)
	id <- substr(range.object$id[[1]], 1, 5)
	panel.grid(h = 10, v = 10, col = "grey", lty = 3)
	panel.xyplot(x, y, col=col[subscripts], fill=fill[subscripts], ...)
	if(highlight){
		st <- start(range.object)/1e6
		en <- end(range.object)/1e6
		lrect(xleft=st,
		      xright=en,
		      ybottom=ylimit[1],
		      ytop=ylimit[2], col=highlight.col, ...)
	}
	if("other.ranges" %in% names(list(...))){
		ranges <- list(...)[["other.ranges"]]
		if(nrow(ranges) > 0){
			st <- start(ranges)/1e6
			en <- end(ranges)/1e6
			lrect(xleft=st,
			      xright=en,
			      ybottom=ylimit[1],
			      ytop=ylimit[2], col=highlight.col, ...)
		}
	}
	##panel.abline(h=0)
	if(add.segments){
		what <- switch(paste("p", panel.number(), sep=""),
			       p1="distance",
			       p2="offspring",
			       p3="mother",
			       p4="father",
			       NULL)
		stopifnot(!is.null(what))
		if(what=="father")
			cbs.sub <- cbs.segs[cbs.segs$id==paste(id, "03", sep="_"), ]
		if(what=="mother")
			cbs.sub <- cbs.segs[cbs.segs$id==paste(id, "02", sep="_"), ]
		if(what=="offspring")
			cbs.sub <- cbs.segs[cbs.segs$id==paste(id, "01", sep="_"), ]
		if(what=="distance"){
			cbs.sub <- dranges[substr(dranges$id, 1, 5) %in% id, ]
			cbs.sub$seg.mean <- -1*cbs.sub$seg.mean
		}
		if(nrow(cbs.sub) > 0){
			cbs.sub$seg.mean[cbs.sub$seg.mean < ylimit[1]] <- ylimit[1] + 0.2
			cbs.sub$seg.mean[cbs.sub$seg.mean > ylimit[2]] <- ylimit[2] - 0.2
			stopifnot(nrow(cbs.sub) > 0)
			panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2,col="black")#gp=gpar("lwd"=2))
		}
	}
}


amp22fun <- function(pruned, ranges.bf2, trioList, verbose=TRUE,...){
	amp22 <- list()
	for(i in 1:nrow(pruned)){
		df <- constructTrioSetFrom(ranged.data=pruned[i,],
					   trioList=trioList, FRAME=1e6,
					   index.in.chromosome=TRUE,
					   verbose=verbose)
		(amp22[[i]] <- xyplot(logR ~ x | subject,
				      data=df,
				      layout=c(1,4),
				      index.cond=list(4:1),
				      pch=21,
				      cex=0.4,
				      xlab="Mb",
				      ylab="",
				      fill="lightblue",
				      border="grey60",
				      range.object=pruned[i,],
				      dranges=ranges.bf2[ranges.bf2$id==pruned$id[i], ],
				      panel=logrpanelfunction2,
				      alpha=1,
				      highlight=FALSE, main=pruned$id[i], ...))
	}
	return(amp22)
}

gridplot <- function(rd, trioList, ylim=c(-2.5,1),
		     highlight=FALSE,
		     alpha=1,
		     highlight.fill="lightblue",
		     highlight.col="grey30",
		     highlight.border="grey60",
		     highlight.alpha=0.4,
		     pch.col="grey60",
		     cex.genes=0.5,
		     pch=21,
		     cex=0.4,
		     cex.scale=0.5,
		     FRAME=1e6, unit="Mb", ...){
	f1 <- list()
	chr.name <- paste("chr", rd$chrom[[1]], sep="")
	for(i in 1:nrow(rd)){
		df <- constructTrioSetFrom(ranged.data=rd[i,],
					   trioList=trioList,
					   FRAME=FRAME,
					   index.in.chromosome=TRUE,
					   unit=unit,
					   ...)
		fill <- rep("white", nrow(df))
		ii <- which(df$x*1e6 >= start(rd)[i] & df$x*1e6 <= end(rd)[i])
		fill[ii] <- highlight.fill
		col <- rep(pch.col, nrow(df))
		col[ii] <- highlight.col
		if(any(df$logR < ylim[1],na.rm=TRUE)){
			replace.index <- which(df$logR < ylim[1])
			df$logR[replace.index] <- ylim[1] + 0.2 + runif(length(replace.index), -0.2, 0.2)
		}
		if(any(df$logR > ylim[2], na.rm=TRUE)){
			replace.index <- which(df$logR > ylim[2])
			df$logR[replace.index] <- ylim[2] - 0.2 + runif(length(replace.index), -0.2, 0.2)
		}
		f1[[i]] <- list()
		f1[[i]][[1]] <- xyplot(logR ~ x | subject, data=df,
				       panel=logrpanelfunction2,
				       ylimit=ylim,
				       layout=c(1,4),
					index.cond=list(4:1),
					pch=pch,
					cex=cex,
					xlab="Mb",
					ylab="",
				       ##fill="lightblue",
					border="grey60",
					range.object=rd[i,],
					##this needs to be all of the min-distance ranges
					##dranges=ranges.bf2[ranges.bf2$id==rd$id[i], ],
				       scales=list(x=list(tick.number=10, cex=cex.scale, tck=c(1,0)),
				       alternating=rep(1, 4),
				       y=list(cex=cex.scale, tck=c(1,0))),
				       par.strip.text=list(lines=0.8, cex=0.7),
				       alpha=alpha,
				       highlight=highlight,
				       col=col,
				       fill=fill, ...)##, main=rd$id[i]))
		xlimit <- f1[[i]][[1]]$x.limits
		df$labels <- as.character(df$subject)
		df$labels[df$labels=="distance"] <- "genes"
		df$labels[nrow(df)] <- "CNV"
		df$labels <- factor(df$labels, levels=c("Father", "Mother", "Offspring", "genes", "CNV"), ordered=TRUE)
		data(rf)
		data(cnv)
		rf <- rf[!duplicated(rf$geneName), ]
		rf.chr <- rf[rf$txStart/1e6 <= xlimit[2] & rf$txEnd/1e6 >= xlimit[1] & rf$chrom==chr.name, ]
		cnv.chr <- cnv[cnv$txStart/1e6 <= xlimit[2] & cnv$txEnd/1e6 >= xlimit[1] & cnv$chrom==chr.name, ]
		flatBed <- flatten.bed(rf.chr)
		flatBed$start <- flatBed$start/1e3
		flatBed$stop <- flatBed$stop/1e3
		flatBed.cnv <- flatten.bed(cnv.chr)
		flatBed.cnv$start <- flatBed.cnv$start/1e3
		flatBed.cnv$stop <- flatBed.cnv$stop/1e3
		bed.objects <- list(genes=flatBed, cnv=flatBed.cnv)
		f1[[i]][[2]] <- xyplot(baf ~ x | labels,
				       data=df,
				       pch=pch,
				       cex=cex,
				       layout=c(1,5),
				       index.cond=list(5:1),
				       xlab="",
				       ylab="",
				       panel=mypanelfunction,
				       highlight=highlight,
				       highlight.col=highlight.col,
				       highlight.border=highlight.border,
				       highlight.alpha=highlight.alpha,
				       ranges=rd,
				       bed.objects=bed.objects,
				       rows=5,
				       cex.genes=cex.genes,
				       par.strip.text=list(lines=0.8, cex=cex.scale),
				       scales=list(x=list(tick.number=10, cex=cex.scale,
						   tck=c(1,0), alternating=1),
				                   y=list(cex=cex.scale, tck=c(0,1),
						   alternating=c(0,0,2,2,2))),
				       ylim=c(-0.02,1.02),
				       col=col,
				       fill=fill, ...)
	}
	return(f1)
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


mypanelfunction <- function(x, y,
			    col, fill,
			    highlight,
			    highlight.col,
			    highlight.alpha,
			    highlight.border,
			    alpha,
			    ranges,
			    ylim,
			    cex.genes=0.6,
			    rows,
			    bed.objects,
			    ..., subscripts){
	if(panel.number() > 2){
		panel.grid(h = 10, v = 10, col = "grey", lty = 3)
		panel.xyplot(x, y, col=col[subscripts], fill=fill[subscripts], ...)
	}
	if(panel.number() == 2){
		flatBed <- bed.objects[["genes"]]
		panel.flatbed(flat=flatBed,
			      showIso=FALSE, rows=rows,
			      cex=cex.genes)
	}
	if(panel.number() == 1){
		flatBed <- bed.objects[["cnv"]]
		panel.flatbed(flat=flatBed,
			      showIso=FALSE, rows=rows,
			      cex=cex.genes,
			      col="red")
	}
	if(highlight){
		##trellis.focus("panel", 1, panel.number())
		ylim=c(-10, 10)
		st <- start(ranges)/1e6
		en <- end(ranges)/1e6
		lrect(xleft=st,
		      xright=en,
		      ybottom=ylim[1],
		      ytop=ylim[2],
		      col=highlight.col, ...)
##		      border=highlight.border,
##		      col=highlight.fill,
##		      alpha=highlight.alpha)
	}
	if("other.ranges" %in% names(list(...))){
		ranges <- list(...)[["other.ranges"]]
		if(nrow(ranges) > 0){
			ylim=c(0, 1)
			st <- start(ranges)/1e6
			en <- end(ranges)/1e6
			lrect(xleft=st,
			      xright=en,
			      ybottom=ylim[1],
			      ytop=ylim[2], col=highlight.col, ...)
			##		      border=highlight.border,
			##		      col=highlight.fill,
			##		      alpha=highlight.alpha)
		}
	}
}

data.frame.for.rectangles <- function(ranges.all, palette,verbose=TRUE){
	stopifnot("num.mark" %in% colnames(ranges.all))
	stopifnot("triostate" %in% colnames(ranges.all))
	## important to order by chromosome before splitting
	ranges.all <- ranges.all[order(ranges.all$chrom), ]
	##del.states <- c("332", "331", "321", "231", "431", "341", "432", "342", "441", "442", "421")
	##amp.states <- c("334", "224", "114", "124", "214", "324", "234", "124", "214", "314", "134")
	del.states <- deletionStates()
	amp.states <- duplicationStates()
	##is.denovo <- ranges.all$triostate %in% c(del.states, amp.states)
	is.denovo <- isDenovo(ranges.all$triostate)
	if(!any(is.denovo)) return(NULL)
	if(!all(is.denovo)){
		if(verbose) message("states non-denovo states detected.  Excluding non-denovo states")
		ranges.all <- ranges.all[is.denovo, ]
	}
	hemizygous.states <- offspring.hemizygous()
	homozygous.states <- offspring.homozygous()
	isAmp <- substr(ranges.all$triostate, 3, 3) %in% c("4", "5")
	cols <- rep(NA, nrow(ranges.all))
	cols[ranges.all$triostate %in% homozygous.states] <- palette[1]
	cols[ranges.all$triostate %in% hemizygous.states] <- palette[2]
	cols[isAmp] <- palette[3]
	## height of rectangle
	h <- 0.75
	meanSegment <- apply(cbind(start(ranges.all), end(ranges.all)), 1, mean)
	data(chromosomeAnnotation)
	chr.size <- chromosomeAnnotation[1:22, "chromosomeSize"]
	chrom <- ranges.all$chrom
	chr.size <- chr.size[chrom]
	y <- split(ranges.all$id, chrom)
	y <- lapply(y, function(x){
		tmp <- as.numeric(as.factor(x))
		names(tmp) <- as.character(x)
		tmp
	})
	y <- unlist(y)
	nms2 <- paste(ranges.all$chrom, ranges.all$id, sep=".")
	if(!identical(names(y), nms2)){
		y <- y[match(nms2, names(y))]
		stopifnot(identical(names(y), nms2))
	}
	dat <- data.frame(x0=start(ranges.all)/1e6,
			  x1=end(ranges.all)/1e6,
			  y0=y-h/2,
			  y1=y+h/2,
			  chr=ranges.all$chrom,
			  coverage=ranges.all$num.mark,
			  midpoint=meanSegment/1e6,
			  id=ranges.all$id,
			  chr.size=chr.size/1e6,
			  col=cols,
			  y=y,
			  stringsAsFactors=FALSE)
	dat$chr <- as.factor(dat$chr)
	return(dat)
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

##fileExt <- function(x) substr(x, 7, 8)


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
	##is.father2 <- substr(ids, 7, 8) == "03"
	##is.mother2 <- substr(ids, 7, 8) == "02"
	##is.offspring <- substr(ids, 7,8) == "01"
	## verify length of family is 3
##	fs <- split(is.father, f)
##	ms <- split(is.mother, f)
##	os <- split(is.offspring, f)
##	nf <- sapply(fs, sum)
##	nm <- sapply(ms, sum)
##	no <- sapply(os, sum)
##	is.trio <- nf == 1 & nm == 1 & no == 1
##	complete.trios <- names(fs)[is.trio]
##}

getDups <- function(){
	dupfns <- list.files(system.file("extdata", package="Beaty"), pattern="duplicates", full.names=TRUE)
	cidr.fn <- dupfns[grep("cidr", dupfns)]
	blind.fn <- dupfns[grep("blind", dupfns)]
	dups.cidr <- read.csv(cidr.fn, as.is=TRUE)
	dups.blind <- read.csv(blind.fn, as.is=TRUE)
##	name1 <- dups.blind$Rep1_DNA_Name
##	name2 <- dups.blind$Rep2_DNA_Name
##	sns1 <- substr(name1, 18, nchar(name1))
##	sns2 <- substr(name2, 18, nchar(name2))
##	sns.bsset <- sapply(sampleNames(bsSet), function(x) strsplit(x, ".txt")[[1]][[1]])
##	index1 <- match(sns1, sns.bsset)
##	index2 <- match(sns2, sns.bsset)
	##complete.trios <- sampleNames(bsSet)[!bsSet$is.wga]
	list(cidr=dups.cidr,
	     blind=dups.blind)
}

isWgaFamily <- function(ids, dna) unique(ids[dna=="WGA"])

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

pennfig8 <- function(penn.denovo, CHR=8, pass.qc=TRUE, ylab="", main=""){
	colnames(penn.denovo)[grep("nmarkers", colnames(penn.denovo))] <- "num.mark"
	index <- which(penn.denovo$chrom==CHR & penn.denovo$MAD < 0.3 & penn.denovo$dna != "WGA")
	if(pass.qc){
		index <- which(penn.denovo$chrom==CHR & penn.denovo$pass.qc & penn.denovo$is.denovo)
	} else index <- which(penn.denovo$chrom==CHR & penn.denovo$is.denovo)
	penn.df <- data.frame(id=penn.denovo$id[index],
			      state=penn.denovo$state[index],
			      start=start(penn.denovo)[index]/1e6,
			      end=end(penn.denovo)[index]/1e6,
			      n=penn.denovo$num.mark[index])
	penn.df$id <- factor(penn.df$id, levels=unique(penn.denovo$id[index]), ordered=TRUE)
	penn.df$id <- as.integer(penn.df$id)
	##penn.df$y <- as.integer(penn.df$id)
	palette <- brewer.pal(9, "Set1")[1:3]
	mykey <- simpleKey(c("homo-del", "hemi-del", "amp"), points=FALSE,
			   rectangles=TRUE, col=palette[1:3], space="top")
	pennfig <- xyplot(id~start, penn.df,
			   panel=function(x, y, start, end, coverage, label=FALSE, ..., subscripts){
				   panel.grid(h=-1, v=-1)
				   panel.xyplot(x, y, ..., subscripts)
				   h <- 0.75
				   lrect(xleft=start[subscripts], xright=end[subscripts],
					 ybottom=as.integer(y)-h/2,
					 ytop=as.integer(y)+h/2,...)
				   if(label) ltext(x,y,labels=coverage[subscripts])
			   },
			   scales=list(x=list(tick.number=10)),
			   border=rep("black", nrow(penn.df)),
			   fill=rep("white", nrow(penn.df)),
			   start=penn.df$start,
			   end=penn.df$end,
			   coverage=penn.df$nm,
			   xlab="Mb",
			   ylab=ylab,
			  ylim=range(penn.df$id),
			   par.strip.text=list(lines=0.7, cex=0.6), main=main)
	pennfig
}



getAllBayesFactorRanges <- function(path=beadstudiodir()){
	ranges.md <- list()
	for(i in 1:22){
		fname <- file.path(path, paste("ranges.bf", i, ".rda", sep=""))
		load(fname)
		ranges.md[[i]] <- get("ranges.bf")
		rm(ranges.bf)
	}
	ii <- which(sapply(ranges.md, is, "list"))
	if(length(ii) > 0){
		for(k in seq_along(ii)){
			j <- ii[k]
			tmp <- RangedDataList(ranges.md[[j]])
			tmp2 <- stack(tmp)
			## for some reason, the above operation adds an extra column
			tmp2 <- tmp2[, -ncol(tmp2)]
			ranges.md[[j]] <- tmp2
			rm(tmp, tmp2); gc()
		}
	}
	stopifnot(all(sapply(ranges.md, is, "RangedData")))
	rdlist <- RangedDataList(ranges.md)
	rd <- stack(rdlist)
	rd$state <- argmax2state(rd$argmax)
	return(rd)
}

calculateDenovoFrequency <- function(ranges.md, penn.offspring, bychrom=FALSE){
	## frequency of denovo event
	trio.states <- trioStates()
	##tmp <- trio.states[ranges.md$argmax, ]
	##calls <- paste(paste(tmp[, 1], tmp[, 2], sep=""), tmp[,3], sep="")
	##calls <- gsub("5", "4", calls)
	ranges.md$state <- gsub("5", "4", ranges.md$state)
	##not.normal <- ranges.md$state != "333"
	not.normal <- isDenovo(ranges.md$state)
##	if(bychrom){
##		df <- list()
##		for(i in 1:22) {
##			ii <- which(ranges.md$chrom==i)
##			tmp <- ranges.md[ii, ]
##			nn <- not.normal[ii]
##			df[[i]] <- table(tmp$state[nn])
##		}
##		state <- sapply(df, names)
##		chrom <- rep(1:22, sapply(df, length))
##		df <- data.frame(freq=unlist(df),chrom=chrom, state=unlist(state))
##	} else {
	df <- as.data.frame(table(ranges.md$state[not.normal]))
	##}
	df$method <- "mindist"
	not.normal <- penn.offspring$pass.qc & penn.offspring$nmarkers >= 10 & isDenovo(penn.offspring$state)
##	if(bychrom){
##		df2 <- list()
##		for(i in 1:22) {
##			ii <- which(penn.offspring$chrom==i)
##			tmp <- penn.offspring[ii, ]
##			nn <- not.normal[ii]
##			df2[[i]] <- table(tmp$state[nn])
##		}
##		state <- sapply(df2, names)
##		chrom <- rep(1:22, sapply(df2, length))
##		df2 <- data.frame(freq=unlist(df2),chrom=chrom, state=unlist(state))
##	} else
	df2 <- as.data.frame(table(penn.offspring$state[not.normal]))
	df2$method <- "PennCNV"
	if(!bychrom){
		colnames(df2) <- colnames(df) <- c("state", "freq", "method")
	}
	df <- rbind(df, df2)
	##df$is.denovo <- df$call %in% c(deletionStates(), duplicationStates())
	##df$is.denovo <- isDenovo(df$state)
	##df <- df[isDenovo(df$call), ]
	df$col <- rep("white", nrow(df))
	##df$col[df$is.denovo] <- "blue"
	i1 <- df$method=="mindist" ##& df$is.denovo
	f1 <- df$freq[i1]
	names(f1) <- as.character(df$state[i1])
	i2 <- df$method=="PennCNV"## & df$is.denovo
	f2 <- df$freq[i2]
	names(f2) <- as.character(df$state[i2])
	f1 <- f1[order(names(f1))]
	f2 <- f2[order(names(f2))]
	f1 <- f1[names(f1) %in% names(f2)]
	f2 <- f2[names(f2) %in% names(f1)]
	stopifnot(all.equal(names(f1), names(f2)))
	f1 <- f1[match(names(f1), names(f2))]
	df <- data.frame(mindist=f1, penn=f2)
	rownames(df) <- names(f1)
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

discordant332Figs <- function(file, triolist, index, pch=21, cex=1, ylim=c(-3,1.5),
			      alpha=1,
			      highlight.fill="lightblue",
			      highlight.col="grey30",
			      highlight.border="grey60",
			      pch.col="grey", ...){
	if(!missing(file)) trellis.device(device="pdf", file=file, onefile=FALSE)
	for(i in seq_along(index)){
		j <- index[i]
		penn[j, ]
		##md[[j]]
		CHR <- penn$chrom[j]
		load(file.path(beadstudiodir(), paste("ranges.bf", CHR, ".rda", sep="")))
		##rd <- ranges.bf[ss(ranges.bf$id)=="12307", ]
		##trace(gridplot, browser)
		fig <- gridplot(rd=penn[j, ], trioList=trioList, pch=pch,
				ylim=ylim, cex=cex,
				highlight.fill=highlight.fill,
				highlight.border=highlight.border,
				highlight.col=highlight.col,
				pch.col=pch.col)
		gridsetup(lattice.object=fig, rd=penn[j, ], width=8, height=5)
	}
	if(!missing(file)) dev.off()
	TRUE
}

readMultiplexData <- function(filenames){
	fn <- filenames
	dfList <- vector("list", length(fn))
	##options(warn=2)
	for(i in seq_along(fn)){
		##flds <- count.fields(fn[i], sep="\t")
		##rows <- which(flds==27)
		dat <- read.table(fn[i], header=TRUE, as.is=TRUE, sep="\t", fill=TRUE)
		cols <- c(3:5, 12)
		tmp <- as.matrix(dat[cols])
		tmp <- suppressWarnings(matrix(as.numeric(tmp), nrow(tmp), ncol(tmp)))
		colnames(tmp) <- colnames(dat)[cols]
		probename <- dat$Sample.Name
		df <- data.frame(tmp)
		if(i > 1){
			df$probeid <- sapply(dat$Sample.Name, function(x) strsplit(x, "\\.")[[1]][1])
			df$sampleid <- sapply(dat$Sample.Name, function(x) strsplit(x, "\\.")[[1]][2])
			##stopifnot(df$sampleid[1] == strsplit(basename(fn)[i], "\\.")[[1]][[1]])
		} else {
			df$probeid <- dat$Sample.Name
			df$sampleid <- strsplit(basename(fn)[i], "\\.")[[1]][[1]]
		}
		dfList[[i]] <- df
	}
	df <- do.call("rbind", dfList)
	colnames(df)[5] <- "sampleid"
	colnames(df)[6] <- "probeid"
	return(df)
}

probeCoordinates <- function(){
	##extdata/Taqman pre_designed CNV assay*.xlsx
	fd <- matrix(c(16040192, 16080463,
		       17280294, 17303840,
		       18124226, 18151112,
		       18124226, 18151112,
		       19047911, 19061542,
		       19047911, 19061542,
		       19601714, 19637890,
		       20641403, 20667147), byrow=TRUE, ncol=2)
	fd <- as.data.frame(fd)
	fd$probeids <- c("HS01502562",
			 "HS00122839",
			 "HS02053540",
			 "HS02374030",
			 "HS04085807",
			 "HS04082205",
			 "HS02102421",
			 "HS01512766")
	colnames(fd)[1:2] <- c("st", "en")
	return(fd)
}


mapTaqmanSampleIds <- function(multiplex){
	sample.hash <- matrix(c("19225_01", "TW227-01",
				"19225_02", "TW227-02",
				"19225_03", "TW227-03",
				"22169_01", "wc079-01",
				"22169_02", "wc079-02",
				"22169_03", "wc079-03",
				"22040_01", "wc134-01",
				"22040_02", "wc134-02",
				"22040_03", "wc134-03",
				"21153_01", "pu262-01",
				"21153_02", "pu262-02",
				"21153_03", "pu262-03",
				"21184_01", "pu293-01",
				"21184_02", "pu293-02",
				"21184_03", "pu293-03",
				"19053_01", "tw054-01",
				"19053_02", "tw054-02",
				"19053_03", "tw054-03",
				"21165_01", "pu274-01",
				"21165_02", "pu274-02",
				"21165_03", "pu274-03",
				"21155_01", "pu264-01",
				"21155_02", "pu264-02",
				"21155_03", "pu264-03")
			      , ncol=2, byrow=TRUE)
	colnames(sample.hash) <- c("gwasid", "taqmanid")
	sample.hash <- as.data.frame(sample.hash)
	index <- match(multiplex$sampleid, sample.hash$taqman)
	not.missing.index <- which(!is.na(index))
	index <- index[!is.na(index)]
	gwas.id <- rep(NA, nrow(multiplex))
	gwas.id[not.missing.index] <- as.character(sample.hash$gwasid[index])
	gwas.id
}

updateWithGwasFindings <- function(multiplex, fd){
	multiplex$seg.mean <- NA
	multiplex$triostate <- NA
	query <- IRanges(fd$st, fd$en)
	CHR <- 22
	cbs.segs <- loadRangesCbs(beadstudiodir(), pattern=paste("cbs_chr", CHR, ".rda", sep=""), CHR=CHR, name="cbs.segs")
	fname <- file.path(beadstudiodir(), paste("ranges.bf", CHR, ".rda", sep=""))
	load(fname)
	ranges.bf <- get("ranges.bf")
	## unique defined in mybase package
	ids <- unique(multiplex$gwas.id, na.rm=TRUE)
	for(i in seq_along(ids)){
		rd2 <- cbs.segs[cbs.segs$id==ids[i], ]
		subj <- IRanges(start(rd2), end(rd2))
		mm <- matchMatrix(findOverlaps(query, subj))
		frame <- 50e3
		while(length(unique(mm[, 1])) != length(query)){
			##if(nrow(mm) != length(query)){
			message("There's a gap in the WGA data that doesn't cover one of the markers")
			missing.range <- seq(length=length(query))[-mm[, "query"]]
			## make the missing.range slightly larger
			query2 <- query
			start(query2)[missing.range] <- start(query)[missing.range]-frame
			end(query2)[missing.range] <- end(query)[missing.range]+frame
			mm <- matchMatrix(findOverlaps(query2, subj))
			frame <- frame+50e3
		}
		segmeans2 <- sapply(split(rd2$seg.mean[mm[, 2]], mm[, 1]), mean)
		## each row corresponds to the segment mean for a specific probe for that sample
		## -- put the segment means in the corresponding rows of the table
		##table.index <- which(gwas.df$gwas.id == ids[i])
		table.index <- which(multiplex$gwas.id==ids[i])
		## among the rows given by table.index, which ones correspond to the probes in the query
		row.index <- match(fd$probeid, multiplex$probeid[table.index])
		##row.index is an index of an index
		jj <- table.index[row.index]
		##check that its right
		stopifnot(all(multiplex$gwas.id[jj] == ids[i]))
		stopifnot(all(multiplex$probeid[jj] == fd$probeid))
		multiplex$seg.mean[jj] <- segmeans2

		rd3 <- ranges.bf[ss(ranges.bf$id) == ss(ids[i]), ]
		subj <- IRanges(start(rd3), end(rd3))
		mm <- matchMatrix(findOverlaps(query, subj))
		frame <- 50e3
		##while(length(unique(mm[, 1])) != length(query)){
		while(length(unique(mm[, 1])) != length(query)){
			missing.range <- seq(length=length(query))[-mm[, "query"]]
			query2 <- query
			start(query2)[missing.range] <- start(query)[missing.range]-frame
			end(query2)[missing.range] <- end(query)[missing.range]+frame
			mm <- matchMatrix(findOverlaps(query2, subj))
			frame <- frame+50e3
		}
		argmax <- sapply(split(rd3$argmax[mm[,2]], mm[, 1]), function(x) paste(unique(x), collapse="-"))
		tmp <- trioStateNames()
		multiplex$triostate[jj] <- tmp[as.integer(argmax)]
	}
	multiplex$family <- ss(multiplex$gwas.id)
	return(multiplex)
}

taqmanFigs <- function(multiplex, xlim=c(0,5), ylim=c(0,5)){
	tmp <- multiplex[!is.na(multiplex$family) & !is.na(multiplex$triostate), ]
	df <- tmp[tmp$triostate==332, ]
	df2 <- tmp[tmp$triostate==333, ]
	taqmanfigs <- list()
	index1 <- which(multiplex$CN.Calculated > ylim[2])
	index2 <- which(multiplex$seg.mean > ylim[2])
	if(length(index1) > 0)
		multiplex$CN.Calculated[index1] <- ylim[2]
	if(length(index2) > 0)
		multiplex$seg.mean[index2] <- ylim[2]
	taqmanfigs[[1]] <- xyplot(CN.Calculated~seg.mean | probeid + family,
				  df,
				  panel=function(x,y, probeid, sampleid, familyid, ..., subscripts){
					  panel.grid(h=5,v=5)
					  panel.abline(0,1, col="grey", lty=2)
					  col <- makeTransparent("grey")
					  panel.xyplot(x, y, pch=21, cex=1.5, col=col, ...)
					  subs <- function(x) substr(x, 7, 8)
					  ids <- sampleid[subscripts]
					  F <- which(subs(ids) == "03")
					  M <- which(subs(ids) == "02")
					  O <- which(subs(ids) == "01")
					  if(length(F) > 0){
						  panel.points(x[F], y[F], pch="F", col="black", cex=1)
					  }
					  if(length(M) > 0){
						  panel.points(x[M], y[M], pch="M", col="black", cex=1)
					  }
					  if(length(O) > 0){
						  panel.points(x[O], y[O], pch="O", col="blue", cex=1)
					  }
					  sampleid[subscripts]
				  }, probeid=df$probeid, sampleid=df$sampleid,
				  familyid=df$family,
				  par.strip.text=list(cex=0.6),
				  xlab="segment mean from CBS",
				  ylab="TaqMan estimate",
				  ##main="Candidate de novo hemizygous deletions",
				  xlim=xlim,
				  aspect="xy",
				  ylim=ylim)
	taqmanfigs[[2]] <- xyplot(CN.Calculated~seg.mean | probeid + family,
				  df2,
				  panel=function(x,y, probeid, sampleid,
				  familyid, ..., subscripts){
					  ##print(subscripts)
					  panel.grid(h=5,v=5)
					  panel.xyplot(x, y, pch=21, cex=1.2, col="grey", ...)
					  subs <- function(x) substr(x, 7, 8)
					  ids <- sampleid[subscripts]
					  F <- which(subs(ids) == "03")
					  M <- which(subs(ids) == "02")
					  O <- which(subs(ids) == "01")
					  if(length(F) > 0){
						  panel.points(x[F], y[F], pch="F", col="black", cex=0.8)
					  }
					  if(length(M) > 0){
						  panel.points(x[M], y[M], pch="M", col="black", cex=0.8)
					  }
					  if(length(O) > 0){
						  panel.points(x[O], y[O], pch="O", col="black", cex=0.8)
					  }
					  sampleid[subscripts]
				  }, probeid=df2$probeid, sampleid=df2$sampleid,
				  familyid=df2$family,
				  par.strip.text=list(cex=0.6),
				  xlab="segment mean from CBS",
				  ylab="TaqMan estimate", ##main="Normal regions predicted from high-throughput platform")
				  xlim=xlim,
				  aspect="xy",
				  ylim=ylim)
	return(taqmanfigs)
}
load12307 <- function(path=beadstudiodir()){
	load(file.path(path, "ranges.bf8.rda"))
	rd <- RangedDataList(ranges.bf)
	rd <- stack(rd)
	rd <- rd[ss(rd$id)=="12307", ]
	return(rd)
}

printMadFig <- function(madfig, dna){
	pars <- trellis.par.get()
	pars$axis.text$cex <- 0.6
	pars$xlab.text$cex <- 0.8
	pars$ylab.text$cex <- 0.8
	trellis.par.set("axis.text", pars$axis.text)
	trellis.par.set("xlab.text", pars$xlab.text)
	trellis.par.set("ylab.text", pars$xlab.text)
	print(madfig)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	grid.text(x=unit(seq(length=length(table(dna))), "native"), y=unit(0.03, "npc"),
		  label=table(dna), gp=gpar(cex=0.8))#paste("n = ", ns, sep=""))
}

pennDenovoRateBwplot <- function(penn.df){
	tab <- table(penn.df$dna)
	drop.category <- names(tab)[tab  < 10]
	if(length(drop.category) > 0) {
		penn.df <- penn.df[!penn.df$dna %in% drop.category, ]
		tab <- table(penn.df$dna)
	}
	fill <- rep("grey90", length(names(tab)))
	##fill[names(tab) == "HapMap"] <- "grey40"
	fig <- bwplot(log10(number) ~ dna, data=penn.df,
				    panel=function(...){
					    panel.grid(h=10,v=0)
					    panel.bwplot(...)
				    },
				    xlab="DNA source", ylab=expression(log[10]("number denovo")),
				    fill=fill, main = "PennCNV (chr 1-22)")
	return(list(fig=fig,tab=tab))
}

printDenovoRateFig <- function(drfig){
	pars <- trellis.par.get()
	pars$axis.text$cex <- 0.6
	pars$xlab.text$cex <- 0.8
	pars$ylab.text$cex <- 0.8
	pars$main.text$cex <- 0.9
	trellis.par.set("axis.text", pars$axis.text)
	trellis.par.set("xlab.text", pars$xlab.text)
	trellis.par.set("ylab.text", pars$ylab.text)
	trellis.par.set("main.text", pars$main.text)
	print(drfig[[1]])
	trellis.focus("panel", 1,1, highlight=FALSE)
	grid.text(x=unit(seq_along(drfig[[2]]), "native"), y=unit(-0.15, "native"), gp=gpar("cex"=0.9),
		  label=as.integer(drfig[[2]]))
	NULL
}

printDenovoFreq <- function(object,...){
	par(las=1)
	black <- makeTransparent("black", alpha=0.3)
	suppressWarnings(graphics:::plot(log10(as.integer(object$mindist)), log10(as.integer(object$penn)), pch=21,
					 cex=3,
					 xlab="min dist", ylab="PennCNV", #pty="s",
					 main=expression(log[10]("frequency")),
					 col=black,...))
	abline(0, 1, col="grey", lty=2)
	text(log10(as.integer(object$mindist)), log10(as.integer(object$penn)), labels=rownames(object), cex=0.7)
	box(col="grey")
}


printchr8fig <- function(pennfig, pennfig.all){
	pars <- trellis.par.get()
	pars$axis.text$cex <- 0.6
	pars$xlab.text$cex <- 0.8
	pars$ylab.text$cex <- 0.8
	pars$main.text$cex <- 0.9
	trellis.par.set("axis.text", pars$axis.text)
	trellis.par.set("xlab.text", pars$xlab.text)
	trellis.par.set("ylab.text", pars$ylab.text)
	trellis.par.set("main.text", pars$main.text)
	grid.newpage()
	lvp <- viewport(x=0, width=unit(0.5, "npc"), just="left", name="lvp")
	pushViewport(lvp)
	pushViewport(dataViewport(xscale=pennfig.all$x.limits,
				  yscale=c(0,1), clip="on"))
	print(pennfig.all, newpage=FALSE, prefix="plot1", more=TRUE)
	upViewport(0)
	print(pennfig, position=c(0.5, 0, 1, 1), more=TRUE, prefix="afterQC")
	upViewport(0)
	grid.text("Chromosome 8", x=unit(0.5, "npc"), y=unit(0.97, "npc"), gp=gpar("cex"=0.9))
}


getPennOffspringRanges <- function(penn.all, bsSet){
	fmo <- substr(penn.all$id, 8, 8)
	penn.offspring <- penn.all[fmo==1, ]

	complete.ids <- unique(sssampleNames(bsSet)[bsSet$complete.trio])
	##penn.offspring <- penn.offspring[ss(penn.offspring$id) %in% complete.ids, ]
	penn.offspring$complete.trio <- ss(penn.offspring$id) %in% complete.ids
	penn.offspring$state <- harmonizeStates(penn.offspring)
	penn.offspring$is.denovo <- penn.offspring$triostate %in% as.character(c(duplicationStatesPenn(), deletionStates()))
	penn.offspring$dna <- bsSet$dna[match(penn.offspring$id, ssampleNames(bsSet))]
	open(bsSet$MAD)
	penn.offspring$MAD <- bsSet$MAD[match(penn.offspring$id, ssampleNames(bsSet))]
	penn.offspring$pass.qc <- bsSet$pass.qc[match(penn.offspring$id, ssampleNames(bsSet))]
	penn.offspring$plate <- bsSet$Sample.Plate[match(penn.offspring$id, ssampleNames(bsSet))]
	wga.families <- unique(sssampleNames(bsSet)[bsSet$is.wga])
	penn.offspring$dna[ss(penn.offspring$is) %in% wga.families] <- "WGA"
	return(penn.offspring)
}


pennDenovoFreq <- function(penn.offspring, bsSet, MIN.COV=10){
	denovo.index <- which(penn.offspring$is.denovo & nMarkers(penn.offspring) >= MIN.COV)
	is.denovo.split <- split(denovo.index, penn.offspring$id[denovo.index])
	denovo.freq <- sapply(is.denovo.split, length)
	dna.source <- bsSet$dna[match(names(is.denovo.split), ssampleNames(bsSet))]
	plate <- bsSet$Sample.Plate[match(names(is.denovo.split), ssampleNames(bsSet))]
	wga.families <- unique(sssampleNames(bsSet)[bsSet$is.wga])
	dna.source[ss(names(is.denovo.split)) %in% wga.families] <- "WGA"
	open(bsSet$MAD)
	mad.offspring <- bsSet$MAD[match(names(is.denovo.split), ssampleNames(bsSet))]
	close(bsSet$MAD)
	wga <- ifelse(dna.source=="WGA", "WGA", "other")
	##graphics:::plot(dat2$mad, log10(dat2$number))
	dat2 <- data.frame(list(number=denovo.freq,
				dna=dna.source,
				wga=wga,
				mad=mad.offspring,
				plate=plate,
				id=names(is.denovo.split)), stringsAsFactors=FALSE)##sampleNames(bsSet)[index]))
	return(dat2)
}

harmonizeRangedData <- function(penn.object){
	cns <- c("chrom", "nmarkers",  "id",
		 "triostate", "pedId", "filename",
		 "complete.trio", "state", "is.denovo",
		 "dna", "MAD", "pass.qc", "plate")
	stopifnot(all.equal(cns, colnames(penn.object)))
	message("Use num.mark for coverage")
	cns[2] <- "num.mark"
	colnames(penn.object) <- cns
	message("Dropping 'pedId'")
	penn.object <- penn.object[, -match("pedId", colnames(penn.object))]
	##message("Harmonizing triostates")
	##penn.object$state <- harmonizeStates(penn.object)
	penn.object
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

chromosomeFig <- function(penn.rd, md.rd, CHR, state, palette=brewer.pal(3, "Set2"), verbose=TRUE, sampleLabels.cex=0.5, ...){
	md.rd <- md.rd[md.rd$chrom==CHR & md.rd$state != "333", ]
	penn.rd <- penn.rd[penn.rd$chrom==CHR & penn.rd$state != "333", ]
	penn.rd <- penn.rd[penn.rd$num.mark >= 10, ]
	penn.rd$triostate <- penn.rd$state
	penn.rd <- penn.rd[isDenovo(penn.rd$state), ]
	md.rd$triostate <- md.rd$state
		if(!missing(state)){
		penn.rd <- penn.rd[penn.rd$state==state,]
		md.rd <- md.rd[md.rd$state==state, ]
	}
	if(nrow(md.rd) > 0){
		md.df <- data.frame.for.rectangles(md.rd, palette=palette, verbose=verbose)
	} else md.df <- NULL
	##penn.denovo <- penn.offspring[penn.offspring$pass.qc & penn.offspring$num.mark >= 10 & penn.offspring$chrom == 22 & isDenovo(penn.offspring$state), ]
	penn.df <- data.frame.for.rectangles(penn.rd, palette, verbose=verbose)
	##penn.df$id <- ss(penn.df$id)
	if(!is.null(md.df)) df <- rbind(md.df, penn.df) else df <- penn.df
	df$id <- ss(df$id)
	df$y <- as.integer(as.factor(df$id))
	if(!is.null(md.df)) {
		df$method <- rep(c("md", "penn"), c(nrow(md.df), nrow(penn.df)))
	} else df$method <- "penn"
	grid.newpage()
	if(length(unique(df$id)) > 45){
		sampleLabels.cex <- 0.4
		if(length(unique(df$id)) > 75)
			sampleLabels.cex <- 0.2
	}
	fig <- p1(df=df, palette=palette, sampleLabels.cex=sampleLabels.cex, ...)
	return(fig)
}

gridwrap <- function(trioList, rd, index, ylim=c(-2.5,1), unit="Mb", ...){
	fig <- gridplot(rd=rd[index, ],
			trioList=trioList, highlight.col="royalblue",
			highlight.border="grey30", pch=".", cex=2,
			ylim=ylim,
			unit=unit,
			...)
	gridsetup(lattice.object=fig, rd=rd[index, ], width=8, height=5)
}
