who <- function(x){
	##exclude CIDR
	y <- rep(NA, length(x))
	y[grep("CIDR", x)] <- "CIDR"
	x <- substr(x, 7, 8)
	x <- paste("id", x, sep="")
	y[x == "id01"] <- "offspring"
	y[x == "id02"] <- "mother"
	y[x == "id03"] <- "father"
	y[x == "id04" | x == "id05" | x == "id06"] <- "?"
	return(y)
}



isInformative <- function(object){
	whoisit <- sapply(object[[1]]$familyMember, who)
	father <- which(whoisit == "father")
	mother <- which(whoisit == "mother")
	offspring <- which(whoisit == "offspring")
	##0. parents homozygous for opposite alleles
	informative.0 <- (calls(object)[, father] == 1 & calls(object)[, mother] == 3) | (calls(object)[, father] == 3 & calls(object)[, mother] == 1)
	#the following are informative only if it is nonbiparental
	##1. F AA, M AB, O !AA
	informative.1 <- calls(object)[, father] == 1 & calls(object)[, mother] == 2 & calls(object)[, offspring] != 1
	##2. M AA, F AB, O !AA
	informative.2 <- calls(object)[, mother] == 1 & calls(object)[, father] == 2 & calls(object)[, offspring] != 1
	##3. F BB, M AB, O !BB
	informative.3 <- calls(object)[, father] == 3 & calls(object)[, mother] == 2 & calls(object)[, offspring] != 3
	##4. M BB, F AB, O !BB
	informative.4 <- calls(object)[, mother] == 3 & calls(object)[, father] == 2 & calls(object)[, offspring] != 3
	informative <- informative.0 | informative.1 | informative.2 | informative.3 | informative.4
	informative[is.na(informative)] <- FALSE
	informative
}

offspringHeterozygous <- function(object, isInf){
	whoisit <- sapply(object[[1]]$familyMember, who)
##	father <- which(whoisit == "father")
##	mother <- which(whoisit == "mother")
	offspring <- which(whoisit == "offspring")
	heterozygous <- calls(object)[, offspring] == 2
	heterozygous[is.na(heterozygous)] <- FALSE
	return(heterozygous)
}

##the probability of observing genotypes consistent with biparental inheritance when in truth the alleles were inherited from only one parent is 0.5
##the probability of observing genotypes consistent with biparental inheritance when in truth the alleles were inherited from both parents is 0.999
##the probability of observing genotypes not consistent with biparental inheritance when in truth the alleles were inherited from both parents is 1-0.999 -- this corresponds to the genotyping error
##(we could potentially have genotyping error estimated from the crlmm confidence score)
##emission.genotypes <- function(object,
##			       prGtError=c(0.001, 0.1), isInf, isHet){##, prBPI, isInf, isHet){
##	if(length(prBPI) != 2) stop("must specific Pr(BPI | hidden state) -- two values.")
##
##	##two states
##	## - by default, states have equal probability
##	e <- matrix(0, nrow(object), 2)
##
##	## hidden state: biparental inheritance
##	e[isInf & !isHet, 1] <- prGtError[1]  ##probability of genotyping error given truth is biparental
##	e[isInf & isHet, 1] <- 1-prGtError[1]
##
##	## hidden state: nonbiparental inheritance
##	e[isInf & !isHet, 2] <- 1-prGtError[2] ## child observed homozygous, true state is biparental
##	e[isInf & isHet, 2] <-  prGtError[2]  ##if true state is BPI, the parents are informative and child is observed heterozygous it must be a genotyping error
##
##	e[!isInf, ] <- 0
####	e[!isInf, 1] <- priorProbBiparental
####	e[!isInf & !isHet, ] <- 0  #child is homozygous.  parental genotypes not informative.  Give equal probability to both states a priori
####	e[!isInf & isHet, 1] <- 1-prGtError[1]  #child is heterozygous and parental genotypes not informative. Give a prior probability of 1-the genotyping error for normal state
####	e[!isInf & isHet, 2] <- prGtError[2]  #child is heterozygous -- must be a genotyping error if inheritance is not biparental
##
##	colnames(e) <- c("Pr(observed | BPI)", "Pr(observed | not BPI)")
##	return(e)
##}
crlmmIlluminaRS2 <- function(sampleSheet=NULL,
			    arrayNames=NULL,
			    batch,
			    ids=NULL,
			    path=".",
			    arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
			    highDensity=FALSE,
			    sep="_",
			    fileExt=list(green="Grn.idat", red="Red.idat"),
			    stripNorm=TRUE,
			    useTarget=TRUE,
			    row.names=TRUE,
			    col.names=TRUE,
			    probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
			    seed=1, save.ab=FALSE, snpFile, cnFile,
			    mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
			    cdfName, sns, recallMin=10, recallRegMin=1000,
			    returnParams=FALSE, badSNP=.7,
			    copynumber=FALSE,
			    load.it=TRUE) {
	if(missing(cdfName)) stop("must specify cdfName")
	if(!isValidCdfName(cdfName)) stop("cdfName not valid.  see validCdfNames")
	if(missing(sns)) sns <- basename(arrayNames)
	if(missing(batch)){
		warning("The batch variable is not specified. The scan date of the array will be used as a surrogate for batch.  The batch variable does not affect the preprocessing or genotyping, but is important for copy number estimation.")
	} else {
		if(length(batch) != length(sns))
			stop("batch variable must be the same length as the filenames")
	}
	batches <- splitIndicesByLength(seq(along=arrayNames), ocSamples())
	k <- 1
	for(j in batches){
		if(verbose) message("Batch ", k, " of ", length(batches))
		RG <- readIdatFiles(sampleSheet=sampleSheet[j, ],
				     arrayNames=arrayNames[j],
				     ids=ids,
				     path=path,
				     arrayInfoColNames=arrayInfoColNames,
				     highDensity=highDensity,
				     sep=sep,
				     fileExt=fileExt,
				     saveDate=TRUE)
		RG <- RGtoXY(RG, chipType=cdfName)
		protocolData <- protocolData(RG)
		res <- preprocessInfinium2(RG,
					   mixtureSampleSize=mixtureSampleSize,
					   fitMixture=TRUE,
					   verbose=verbose,
					   seed=seed,
					   eps=eps,
					   cdfName=cdfName,
					   sns=sns[j],
					   stripNorm=stripNorm,
					   useTarget=useTarget)
		rm(RG); gc()
		## MR: number of rows should be number of SNPs + number of nonpolymorphic markers.
		##  Here, I'm just using the # of rows returned from the above function
		if(k == 1){
			if(verbose) message("Initializing container for alleleA, alleleB, call, callProbability")
			load.obj <- loadObject("callSet", load.it)
			if(!load.obj){
				callSet <- new("SnpSuperSet",
					       alleleA=initializeBigMatrix(name="A", nr=nrow(res[[1]]), nc=length(sns)),
					       alleleB=initializeBigMatrix(name="B", nr=nrow(res[[1]]), nc=length(sns)),
					       call=initializeBigMatrix(name="call", nr=nrow(res[[1]]), nc=length(sns)),
					       callProbability=initializeBigMatrix(name="callPr", nr=nrow(res[[1]]), nc=length(sns)),
					       annotation=cdfName)
				sampleNames(callSet) <- sns
				save(callSet, file=file.path(ldPath(), "callSet.rda"))
			} else load(file.path(ldPath(), "callSet.rda"))
			phenoData(callSet) <- getPhenoData(sampleSheet=sampleSheet,
							   arrayNames=sns,
							   arrayInfoColNames=arrayInfoColNames)
			pD <- data.frame(matrix(NA, length(sns), 1), row.names=sns)
			colnames(pD) <- "ScanDate"
			protocolData(callSet) <- new("AnnotatedDataFrame", data=pD)
			pData(protocolData(callSet))[j, ] <- pData(protocolData)
			featureNames(callSet) <- res[["gns"]]
			pData(callSet)$SNR <- initializeBigVector("crlmmSNR-", length(sns), "double")
			pData(callSet)$SKW <- initializeBigVector("crlmmSKW-", length(sns), "double")
			pData(callSet)$gender <- rep(NA, length(sns))
			mixtureParams <- initializeBigMatrix("crlmmMixt-", nr=4, nc=ncol(callSet), vmode="double")
			save(mixtureParams, file=file.path(ldPath(), "mixtureParams.rda"))
			if(missing(batch)){
				protocolData(callSet)$batch <- rep(NA, length(sns))
			} else{
				protocolData(callSet)$batch <- batch
			}
			featureData(callSet) <- addFeatureAnnotation(callSet)
			open(mixtureParams)
			open(callSet$SNR)
			open(callSet$SKW)
		}
		if(k > 1 & nrow(res[[1]]) != nrow(callSet)){
			##RS: I don't understand why the IDATS for the
			##same platform potentially have different lengths
			res[["A"]] <- res[["A"]][res$gns %in% featureNames(callSet), ]
			res[["B"]] <- res[["B"]][res$gns %in% featureNames(callSet), ]
		}
		if(missing(batch)){
			protocolData(callSet)$batch[j] <- as.numeric(as.factor(protocolData$ScanDate))
		}
		## MR: we need to define a snp.index vs np.index
		snp.index <- match(res$gns, featureNames(callSet))
		A(callSet)[snp.index, j] <- res[["A"]]
		B(callSet)[snp.index, j] <- res[["B"]]
		pData(callSet)$SKW[j] <- res$SKW
		pData(callSet)$SNR[j] <- res$SNR
		mixtureParams[, j] <- res$mixtureParams
		rm(res); gc()
		k <- k+1
	}
	return(callSet)
}


getFamilyId <- function(samplesheet){
	familyId <- sapply(samplesheet[, 4], function(x) strsplit(x, "_")[[1]][[1]])
}


setMethod("addSampleSheet", "SnpSuperSet", function(object){
	data(samplesheet)
	samplesheet <- get("samplesheet")

	pd1 <- pData(object)
	pd2 <- samplesheet[match(sampleNames(object), samplesheet$sampleNames), ]
	pd <- cbind(pd1, pd2)
	pD <- new("AnnotatedDataFrame",
		  data=pd,
		  varMetadata=data.frame(labelDescription=colnames(pd),
		  row.names=colnames(pd)))
	return(pD)
})


addSampleSheetAnnotation <- function(object, samplesheet){
	if("familyId" %in% varLabels(object)) return()
	pD <- pData(object[[1]])
	pd2 <- samplesheet[match(sampleNames(object), samplesheet$sampleNames), ]
	##The trio information is in the sample name:
	##excludeIndex <- grep("CIDR", pd2[, 4])
	##crlmmSetList <- crlmmSetList[, -excludeIndex]

	familyId <- sapply(pd2[, 4], function(x) strsplit(x, "_")[[1]][[1]])
	familyMember <- sapply(pd2[, 4], function(x) strsplit(x, "_")[[1]][[2]])
	##table(familyMember, pd2[, "Gender"])
	## 01: offspring
	## 02: Mother
	## 03: Father
	pd2 <- cbind(pD, pd2, as.character(familyId), as.character(familyMember))
	colnames(pd2)[34:35] <- c("familyId", "familyMember")
	pD <- new("AnnotatedDataFrame", data=pd2, varMetadata=data.frame(labelDescription=colnames(pd2)))
	return(pD)
##	phenoData(crlmmSetList[[1]]) <- pD
	pData(crlmmSetList[[1]])[, "familyId"] <- as.character(pData(crlmmSetList[[1]])[, "familyId"])
	pData(crlmmSetList[[1]])[, "familyMember"] <- as.character(pData(crlmmSetList[[1]])[, "familyMember"])
}

scaleSnr <- function(snr, min.scale, max.scale){
	tmp <- (snr-min(snr))/(max(snr)-min(snr))  ## 0 -> 1
	##max.scale <- 1
	##min.scale <- 0
	b <- 1/(max.scale-min.scale)
	a <- min.scale*b
	bg.scale <- (tmp + a)/b
	return(bg.scale)
}

sourceCrlmm <- function(){
	library(genefilter); library(affyio)
	.crlmmPkgEnv <- new.env()
	source("~/madman/Rpacks/crlmm/R/AllClasses.R")
	source("~/madman/Rpacks/crlmm/R/AllGenerics.R")
	source("~/madman/Rpacks/crlmm/R/cnrma-functions.R")
	source("~/madman/Rpacks/crlmm/R/crlmm-functions.R")
	source("~/madman/Rpacks/crlmm/R/crlmm-illumina.R")
	source("~/madman/Rpacks/crlmm/R/utils.R")
	source("~/madman/Rpacks/crlmm/R/methods-ABset.R")
	source("~/madman/Rpacks/crlmm/R/methods-eSet.R")
	source("~/madman/Rpacks/crlmm/R/methods-CopyNumberSet.R")
	source("~/madman/Rpacks/crlmm/R/methods-CrlmmSet.R")
	source("~/madman/Rpacks/crlmm/R/methods-SnpCallSetPlus.R")
}

pennCnv2matrix <- function(object, CHR, cdfName, platform){
	stopifnot(isValidCdfName(cdfName, platform))
	pkgname <- paste(cdfName, "Crlmm", sep="")
	path <- system.file("extdata", package=pkgname)
	loader("cnProbes.rda", pkgname=pkgname, envir=.crlmmPkgEnv)
	cnProbes <- get("cnProbes", envir=.crlmmPkgEnv)
	loader("snpProbes.rda", pkgname=pkgname, envir=.crlmmPkgEnv)
	snpProbes <- get("snpProbes", envir=.crlmmPkgEnv)
	snpProbes <- snpProbes[snpProbes$chr == CHR, ]
	cnProbes <- cnProbes[cnProbes$chr==CHR, ]
	probes <- rbind(snpProbes, cnProbes)
	probes <- probes[order(probes$position), ]

	state <- object$TrioState
	stateIndex <- NA
	familyMember <- unique(object$FamilyMember)
	if(length(familyMember) > 1) stop("should just be one sample")
	stateIndex <- switch(familyMember,
			     father=1,
			     mother=2,
			     offspring=3,
			     NA)
	if(is.na(stateIndex)) stop("FamilyMember label is not father, mother, or offspring")
	cn.state <- paste("s", substr(state, stateIndex, stateIndex), sep="")
	getCnState <- function(x){
		switch(x,
		       s1=0,
		       s2=1,
		       s3=2,
		       s4=2,
		       s5=3,
		       s6=4,
		       NA)
	}
	copynumber <- sapply(cn.state, getCnState)
	if(any(is.na(copynumber))) warning("one or more of the assigned copy number states is missing")

	len <- object$LengthCNV
	len <- as.integer(sapply(len, function(x) paste(unlist(strsplit(x,",")), collapse="")))
	nsnps <- object$NumberSNPs
	x <- rep(as.integer(2), nrow(probes))
	for(i in 1:nrow(object)){
		cond <- probes$position >= object$StartPosition[i] & probes$position <= object$EndPosition[i]
		x[cond] <- as.integer(copynumber[i])
	}
	names(x) <- rownames(probes)
	return(x)
}








posteriorNonpolymorphic <- function(plateIndex, envir, priorProb, cnStates=0:6){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	if(missing(priorProb)) priorProb <- rep(1/length(cnStates), length(cnStates)) ##uniform
	plate <- envir[["plate"]]
	uplate <- envir[["plate"]]
	NP <- envir[["NP"]][, plate==uplate[p]]
	nuT <- envir[["nuT"]][, p]
	phiT <- envir[["phiT"]][, p]
	sig2T <- envir[["sig2T"]][, p]
	##Assuming background variance for np probes is the same on the log-scale
	emit <- array(NA, dim=c(nrow(NP), ncol(NP), length(cnStates)))##SNPs x sample x 'truth'
	lT <- log2(NP)
	sds <- sqrt(sig2T)
	counter <- 1##state counter
	for(CT in cnStates){
		cat(".")
		if(CHR == 23) browser()
		means <- suppressWarnings(log2(nuT + CT*phiT))
		emit[, , counter] <- dnorm(lT, mean=means, sd=sds)
		counter <- counter+1
	}
	for(j in seq(along=cnStates)){
		emit[, , j] <- priorProb[j]*emit[, , j]
	}
	homDel <- emit[, , 1]
	hemDel <- emit[, , 2]
	norm <- emit[, , 3]
	amp <- emit[, , 4]
	amp4 <- emit[, , 5]
	amp5 <- emit[, , 6]
	amp6 <- emit[, , 7]
	total <- homDel+hemDel+norm+amp+amp4+amp5+amp6
	weights <- array(NA, dim=c(nrow(NP), ncol(NP), length(cnStates)))
	weights[, , 1] <- homDel/total
	weights[, , 2] <- hemDel/total
	weights[, , 3] <- norm/total
	weights[, , 4] <- amp/total
	weights[, , 5] <- amp4/total
	weights[, , 6] <- amp5/total
	weights[, , 7] <- amp6/total
	##posterior mode
	posteriorMode <- apply(weights, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	posteriorMode <- posteriorMode-1
	##sns <- envir[["sns"]]
	##colnames(posteriorMode) <- sns
	##envir[["np.posteriorMode"]] <- posteriorMode
	##envir[["np.weights"]] <- weights
	posteriorMeans <- 0*homDel/total + 1*hemDel/total + 2*norm/total + 3*amp/total + 4*amp4/total + 5*amp5/total + 6*amp6/total
	##colnames(posteriorMeans) <- sns
	##envir[["np.posteriorMeans"]] <- posteriorMeans
	return(posteriorMode)
}

posteriorWrapper <- function(envir){
	snp.PM <- matrix(NA, length(envir[["snps"]]), length(envir[["sns"]]))
	np.PM <- matrix(NA, length(envir[["cnvs"]]), length(envir[["sns"]]))
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	for(p in seq(along=uplate)){
		tmp <- posteriorPolymoprhic(plateIndex=p, envir=envir)
		snp.PM[, plate==uplate[p]] <- tmp
		##snp.pm <- env[["posteriorMode"]]
		##trace(posteriorNonpolymorphic, browser)
		tmp <- posteriorNonpolymorphic(plateIndex=p, envir=envir)
		np.PM[, plate==uplate[p]] <- tmp##env[["np.posteriorMode"]]
		##pMode <- rbind(snp.pm, np.pm)
		##rownames(pMode) <- c(env[["snps"]], env[["cnvs"]])
		##dn <- dimnames(pMode)
		##pMode <- matrix(as.integer(pMode), nrow(pMode), ncol(pMode))
	}
	PM <- rbind(snp.PM, np.PM)
	PM <- matrix(as.integer(PM), nrow(PM), ncol(PM))
	dns <- list(c(envir[["snps"]], envir[["cnvs"]]), envir[["sns"]])
	dimnames(PM) <- dns
	return(PM)
}


##for polymorphic probes
posteriorPolymorphic <- function(plateIndex, envir, priorProb, cnStates=0:6){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	if(missing(priorProb)) priorProb <- rep(1/length(cnStates), length(cnStates)) ##uniform
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	A <- envir[["A"]]
	B <- envir[["B"]]
	A <- A[, plate==uplate[p]]
	B <- B[, plate==uplate[p]]
	calls <- envir[["calls"]]
	calls <- calls[, plate==unique(plate)[p]]
	probA <- sqrt(rowMeans(calls == 1, na.rm=TRUE))
	probB <- sqrt(rowMeans(calls == 3, na.rm=TRUE))
	sig2A <- envir[["sig2A"]]
	sig2B <- envir[["sig2B"]]
	tau2A <- envir[["tau2A"]]
	tau2B <- envir[["tau2B"]]
	corrA.BB <- envir[["corrA.BB"]]
	corrB.AA <- envir[["corrB.AA"]]
	corr <- envir[["corr"]]
	nuA <- envir[["nuA"]]
	nuB <- envir[["nuB"]]
	phiA <- envir[["phiA"]]
	phiB <- envir[["phiB"]]
	emit <- array(NA, dim=c(nrow(A), ncol(A), 28))##SNPs x sample x 'truth'
	##AAAA, AAAB, AABB, ABBB, BBBB
	##AAAAA, AAAAB, AAABB, AABBB, ABBBB, BBBBB
	##AAAAAA, AAAAAB, AAAABB, AAABBB, AABBBB, ABBBBB, BBBBBB
	lA <- log2(A)
	lB <- log2(B)
	X <- cbind(lA, lB)
	counter <- 1##state counter
	for(CT in cnStates){
		cat(".")
		for(CA in 0:CT){
			CB <- CT-CA
			A.scale <- sqrt(tau2A[, p]*(CA==0) + sig2A[, p]*(CA > 0))
			B.scale <- sqrt(tau2B[, p]*(CB==0) + sig2B[, p]*(CB > 0))
			scale <- c(A.scale, B.scale)
			if(CA == 0 & CB == 0) rho <- 0
			if(CA == 0 & CB > 0) rho <- corrA.BB[, p]
			if(CA > 0 & CB == 0) rho <- corrB.AA[, p]
			if(CA > 0 & CB > 0) rho <- corr[, p]
			if(CHR == 23) browser()
			means <- cbind(suppressWarnings(log2(nuA[, p]+CA*phiA[, p])), suppressWarnings(log2(nuB[, p]+CB*phiB[, p])))
			covs <- rho*A.scale*B.scale
			A.scale2 <- A.scale^2
			B.scale2 <- B.scale^2
			##ensure positive definite
			##Sigma <- as.matrix(nearPD(matrix(c(A.scale^2, covs,
			##covs, B.scale^2), 2, 2))[[1]])
			m <- 1##snp counter
			for(i in 1:nrow(A)){
				Sigma <- matrix(c(A.scale2[i], covs[i], covs[i], B.scale2[i]), 2,2)
				xx <- matrix(X[i, ], ncol=2)
				tmp <- dmvnorm(xx, mean=means[i, ], sigma=Sigma)
				##Using HWE: P(CA=ca, CB=cb|CT=c)
				ptmp <- (probA[i]^CA)*(probB[i]^CB)*tmp
				emit[m, , counter] <- ptmp
				m <- m+1
			}
			counter <- counter+1
		}
	}
	##priorProb=P(CT=c)
	homDel <- priorProb[1]*emit[, , 1]
	hemDel <- priorProb[2]*emit[, , c(2, 3)] # + priorProb[3]*emit[, c(4, 5, 6)] + priorProb[4]*emit[, c(7:10)]
	norm <- priorProb[3]*emit[, , 4:6]
	amp <- priorProb[4]*emit[, , 7:10]
	amp4 <- priorProb[5]*emit[, , 11:15]
	amp5 <- priorProb[6]*emit[, , 16:21]
	amp6 <- priorProb[7]*emit[, , 22:28]
	##sum over the different combinations within each copy number state
	hemDel <- apply(hemDel, c(1,2), sum)
	norm <- apply(norm, c(1, 2), sum)
	amp <- apply(amp, c(1,2), sum)
	amp4 <- apply(amp4, c(1,2), sum)
	amp5 <- apply(amp5, c(1,2), sum)
	amp6 <- apply(amp6, c(1,2), sum)
	total <- homDel+hemDel+norm+amp+amp4+amp5+amp6
	weights <- array(NA, dim=c(nrow(homDel), ncol(A), 7))
	weights[, , 1] <- homDel/total
	weights[, , 2] <- hemDel/total
	weights[, , 3] <- norm/total
	weights[, , 4] <- amp/total
	weights[, , 5] <- amp4/total
	weights[, , 6] <- amp5/total
	weights[, , 7] <- amp6/total
	##posterior mode
	posteriorMode <- apply(weights, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	posteriorMode <- posteriorMode-1
	##This is for one plate.  Need to instantiate a much bigger
	##object in the environment

	##envir[["posteriorMode"]] <- posteriorMode
	##weights <- list(homDel/total, hemDel/total, norm/total, amp/total, amp4/total, amp5/total, amp6/total)
	##names(weights) <- c(cnStates)
	##envir[["weights"]] <- weights
	posteriorMeans <- 0*homDel/total + 1*hemDel/total + 2*norm/total + 3*amp/total + 4*amp4/total + 5*amp5/total + 6*amp6/total
	##sns <- envir[["sns"]]
	##colnames(posteriorMeans) <- sns
	##envir[["posteriorMeans"]] <- posteriorMeans
	return(posteriorMode)
}



##par(mfrow=c(3,3), las=1, pty="s", ask=ask, mar=c(2, 2, 2, 2), oma=c(2, 2, 1, 1))
##indices <- split(snpIndex(crlmmSetList), rep(1:length(snpIndex(crlmmSetList)), each=9, length.out=length(snpIndex(crlmmSetList))))
####for(j in seq(along=indices)[1:10]){
##j <- 1
##	cat(j, "\n")
##	k <- 1
##	for(i in indices[[j]]){
##		gt <- calls(crlmmSetList)[i, ]
##		pch <- as.character(gt)
##		cex <- 0.9
##		plot(crlmmSetList[i, ],
##		     pch=pch,
##		     col=pch.col[gt],
##		     cex=cex,
##		     xlim=xlim, ylim=ylim,
##		     type="n")
##		if(plotpoints){
##			for(b in seq(along=unique(batch(crlmmSetList)))){
##				points(crlmmSetList[i, J[[b]]],
##				       pch=pch,
##				       col=colors[b], bg=colors[b], cex=cex,
##				       xlim=xlim, ylim=ylim)
##			}
##		}
##		for(b in seq(along=unique(batch(crlmmSetList)))){
##			ellipse(crlmmSetList[i, J[[b]]], copynumber=2, col=colors[b], lwd=lwd)
##		}
##		##legend("bottomright", bty="n", legend=featureNames(crlmmSetList)[i])
##		if(k == 1) {
##			legend("bottomleft", bty="n", fill=colors, legend=c("CEPH", "Yoruba", "Asian"))
##			mtext("A", 1, outer=TRUE, line=1)
##			mtext("B", 2, outer=TRUE, line=0)
##		}
##		k <- k+1
##	}

##myIdatReader <- function(sampleSheet=NULL, arrayNames, ids=NULL, path=".",
##			 arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
##			 highDensity=FALSE, sep="_", fileExt=list(green="Grn.idat", red="Red.idat"), saveDate=FALSE){
##	narrays = length(arrayNames)
##	##grnfiles = paste(arrayNames, fileExt$green, sep=sep)
##	grnidats <- grnfiles <- arrayNames[grep("_Grn.idat", arrayNames)]
##	##redfiles = paste(arrayNames, fileExt$red, sep=sep)
##	redidats <- redfiles <- arrayNames[grep("_Red.idat", arrayNames)]
##	if(length(grnfiles)==0 || length(redfiles)==0)
##	       stop("Cannot find .idat files")
##       if(length(grnfiles)!=length(redfiles))
##	       stop("Cannot find matching .idat files")
##	if(!all(c(redfiles,grnfiles) %in% dir(path=path)))
##		stop("Missing .idat files: red\n", paste(redfiles[!(redfiles %in% dir(path=path))], sep=" "), "\n green\n",
##		     paste(grnfiles[!(grnfiles %in% dir(path=path))], sep=" "))
##	##grnidats = file.path(path, grnfiles)
##	##redidats = file.path(path, redfiles)
##
##       headerInfo = list(nProbes = rep(NA, narrays),
##                         Barcode = rep(NA, narrays),
##                         ChipType = rep(NA, narrays),
##                         Manifest = rep(NA, narrays), # not sure about this one - sometimes blank
##                         Position = rep(NA, narrays)) # this may also vary a bit
##       dates = list(decode=rep(NA, narrays),
##                    scan=rep(NA, narrays))
##
##	if(!is.null(sampleSheet)) pd = new("AnnotatedDataFrame", data = sampleSheet)
##
##       # read in the data
##       for(i in seq(along=arrayNames)) {
##	       cat("reading", arrayNames[i], "\t")
##	       idsG = idsR = G = R = NULL
##	       cat(paste(sep, fileExt$green, sep=""), "\t")
##	       G = readIDAT(grnidats[i])
##	       idsG = rownames(G$Quants)
##	       headerInfo$nProbes[i] = G$nSNPsRead
##	       headerInfo$Barcode[i] = G$Barcode
##	       headerInfo$ChipType[i] = G$ChipType
##	       headerInfo$Manifest[i] = G$Unknown$MostlyNull
##	       headerInfo$Position[i] = G$Unknowns$MostlyA
##
##	       if(headerInfo$ChipType[i]!=headerInfo$ChipType[1] || headerInfo$Manifest[i]!=headerInfo$Manifest[1]) {
##		       ## || headerInfo$nProbes[i]!=headerInfo$nProbes[1] ## removed this condition as some arrays used the same manifest
##		       ## but differed by a few SNPs for some reason - most of the chip was the same though
##		       ##           stop("Chips are not of all of the same type - please check your data")
##		       warning("Chips are not of the same type.  Skipping ", basename(grnidats[i]), " and ", basename(redidats[i]))
##		       next()
##	       }
##
##	       dates$decode[i] = G$RunInfo[1, 1]
##	       dates$scan[i] = G$RunInfo[2, 1]
##
##	       if(i==1) {
##		       if(is.null(ids) && !is.null(G)){
##			       ids = idsG
##		       } else  stop("Could not find probe IDs")
##		       nprobes = length(ids)
##		       narrays = length(arrayNames)
##
##		       tmpmat = matrix(NA, nprobes, narrays)
##		       rownames(tmpmat) = ids
##		       if(!is.null(sampleSheet)){
##			       colnames(tmpmat) = sampleSheet$Sample_ID
##		       } else colnames(tmpmat) = arrayNames
##		       RG <- new("NChannelSet",
##				 R=tmpmat, G=tmpmat, Rnb=tmpmat, Gnb=tmpmat,
##				 Rse=tmpmat, Gse=tmpmat, annotation=headerInfo$Manifest[1],
##				 phenoData=pd, storage.mode="environment")
##		       rm(tmpmat)
##		       gc()
##	       }
##
##	       if(length(ids)==length(idsG)) {
##		       if(sum(ids==idsG)==nprobes) {
##			       RG@assayData$G[,i] = G$Quants[, "Mean"]
##			       RG@assayData$Gnb[,i] = G$Quants[, "NBeads"]
##			       RG@assayData$Gse[,i] = G$Quants[, "SD"]
##		       }
##	       }
##	       else {
##		       indG = match(ids, idsG)
##		       RG@assayData$G[,i] = G$Quants[indG, "Mean"]
##		       RG@assayData$Gnb[,i] = G$Quants[indG, "NBeads"]
##		       RG@assayData$Gse[,i] = G$Quants[indG, "SD"]
##	       }
##	       rm(G)
##	       gc()
##
##	       cat(paste(sep, fileExt$red, sep=""), "\n")
##	       R = readIDAT(redidats[i])
##	       idsR = rownames(R$Quants)
##
##	       if(length(ids)==length(idsG)) {
##		       if(sum(ids==idsR)==nprobes) {
##			       RG@assayData$R[,i] = R$Quants[ ,"Mean"]
##			       RG@assayData$Rnb[,i] = R$Quants[ ,"NBeads"]
##			       RG@assayData$Rse[,i] = R$Quants[ ,"SD"]
##		       }
##	       }
##	       else {
##		       indR = match(ids, idsR)
##		       RG@assayData$R[,i] = R$Quants[indR, "Mean"]
##		       RG@assayData$Rnb[,i] = R$Quants[indR, "NBeads"]
##		       RG@assayData$Rse[,i] = R$Quants[indR, "SD"]
##	       }
##	       rm(R)
##	       gc()
##       }
##       if(saveDate) {
##	       protocolData(RG)[["ScanDate"]] = dates$scan
##       }
##       storageMode(RG) = "lockedEnvironment"
##       RG
##}

isBiparental.matrix <- function(object, allowHetParent=TRUE){
	if(!all(colnames(object) == c("father", "mother", "offspring"))) stop()
	F <- object[, 1]
	M <- object[, 2]
	O <- object[, 3]
	##M/F AA, F/M BB, O AB
	##isHet <- offspringHeterozygous(object)  ##offspring is heterozygous
	biparental <- rep(NA, nrow(object))
	biparental[F==1 & M == 3 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	biparental[F==3 & M == 1 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	##M/F AA, F/M BB, O AA or BB
	biparental[F==1 & M == 3 & (O == 1 | O == 3)] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	biparental[F==3 & M == 1 & (O == 1 | O == 3)] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	## M/F AA, F/M BB, O AB
	if(allowHetParent) biparental[F == 1 & M == 2 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	if(allowHetParent) biparental[F == 2 & M == 1 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	## F AB, M AA, O BB is not biparental
	## F AA, M AB, O BB is not biparental
	biparental[F == 2 & M == 1 & O == 3] <- FALSE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	biparental[F == 1 & M == 2 & O == 3] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	## M AA, F AB, O AB
	if(allowHetParent) biparental[F == 2 & M == 3 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	if(allowHetParent) biparental[F == 3 & M == 2 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	## F=AB, M=BB, O=AA is NOT biparental
	biparental[F == 2 & M == 3 & O == 1] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	biparental[F == 3 & M == 2 & O == 1] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	return(biparental)
}



isBiparental.SnpCallSetPlus <- function(object, allowHetParent=TRUE){
	## There will be more noninformative calls when allowHetParent is FALSE.
	## When allowHetParent is FALSE,
	##      i.   F=AA, M=AB, O=AB is treated as noninformative (value NA)
	##      ii.  F=AA, M=AB, O=BB is informative and regarded as nonbiparental inheritance
	##  (i.)   often arises from UPD and would not be a mendelian inconsistency (MI) resulting from a deletion in the offspring
	##  (ii.)  Is an MI that could arise as a result of a denovo deletion
	## Regardless of allowHetParent,
	##       F=AA, M=AB, O=AA is treated as noninformative
	##       F=BB, M=AB, O=BB is treated as noninformative
	## values returned:
	##  NA : no information on whether the inheritance was biparental
	##  TRUE:  inheritance is biparental
	##  FALSE: evidence of non-biparental inheritance
	if(length(object$familyMember) < 3) stop("object$familyMember not the right length")
	whoisit <- sapply(object$familyMember, who)
	father <- which(whoisit == "father")
	mother <- which(whoisit == "mother")
	offspring <- which(whoisit == "offspring")
	F <- calls(object[, father])
	M <- calls(object[, mother])
	O <- calls(object[, offspring])

	object <- cbind(F, M, O)
	colnames(object) <- c("father", "mother", "offspring")
	biparental <- isBiparental.matrix(object, allowHetParent=allowHetParent)
	return(biparental)
##	##M/F AA, F/M BB, O AB
##	##isHet <- offspringHeterozygous(object)  ##offspring is heterozygous
##	biparental <- rep(NA, nrow(object))
##	biparental[F==1 & M == 3 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	biparental[F==3 & M == 1 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	##M/F AA, F/M BB, O AA or BB
##	biparental[F==1 & M == 3 & (O == 1 | O == 3)] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	biparental[F==3 & M == 1 & (O == 1 | O == 3)] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	## M/F AA, F/M BB, O AB
##	if(allowHetParent) biparental[F == 1 & M == 2 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	if(allowHetParent) biparental[F == 2 & M == 1 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	## F AB, M AA, O BB is not biparental
##	## F AA, M AB, O BB is not biparental
##	biparental[F == 2 & M == 1 & O == 3] <- FALSE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	biparental[F == 1 & M == 2 & O == 3] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	## M AA, F AB, O AB
##	if(allowHetParent) biparental[F == 2 & M == 3 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	if(allowHetParent) biparental[F == 3 & M == 2 & O == 2] <- TRUE#Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	## F=AB, M=BB, O=AA is NOT biparental
##	biparental[F == 2 & M == 3 & O == 1] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	biparental[F == 3 & M == 2 & O == 1] <- FALSE#Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	return(biparental)
}
## Twelve conditions for biparental inheritance.
hmm.options <- function(states=c("BPI", "notBPI"),
			initialP=c(0.99, 0.01),
			TAUP=1e7,
			prGtError=c(0.001, 0.01),
			verbose=FALSE,
			allowHetParent=FALSE,
			normalIndex=1){
	names(prGtError) <- states
	names(initialP) <- states
	list(states=states,
	     initialP=initialP,
	     TAUP=TAUP,
	     prGtError=prGtError,
	     verbose=verbose,
	     allowHetParent=allowHetParent,
	     normalIndex=1)
}

isFMOtrio <- function(Sample.Name){ ## is Father, Mother, Offspring trio
	familyId <- substr(Sample.Name, 1, 5)
	if(length(unique(familyId)) > 1) stop("must be one family")
	individualId <- substr(Sample.Name, 7, 8)
	whoisit <- sapply(individualId, who)
	##whoisit <- sapply(object$familyMember, who)
	father <- which(whoisit == "father")
	mother <- which(whoisit == "mother")
	offspring <- which(whoisit == "offspring")
	if(length(father) == 1 & length(mother) == 1 & length(offspring) == 1){
		return(TRUE)
	} else {
		return(FALSE)
	}
}

hmm.SnpCallSetPlus <- function(object, hmmOptions){
	require(VanillaICE) || stop("VanillaICE not available")
	TAUP <- hmmOptions[["TAUP"]]
	states <- hmmOptions[["states"]]
	initialP <- hmmOptions[["initialP"]]
	verbose <- hmmOptions[["verbose"]]
	familyId <- substr(object$Sample.Name, 1, 5)
	object$familyId <- familyId
	trios.complete <- split(object$Sample.Name, familyId)
	trios.complete <- trios.complete[sapply(trios.complete, length) >= 3]
	if(length(trios.complete) == 0){
		warning(paste("No complete trios in plate", unique(object$Sample.Plate)))
		return(NULL)
	}
	isFMO <- sapply(trios.complete, isFMOtrio)
	trios.complete <- trios.complete[isFMO]
	Sample.Names <- unlist(trios.complete)
	object <- object[, object$Sample.Name %in% Sample.Names]
	if(ncol(object) < 3){
		return("No complete trios")
	}
	familyId <- object$familyId
	fit <- matrix(NA, nrow(object), ncol=length(unique(familyId)))
	colnames(fit) <- unique(familyId)
	rownames(fit) <- featureNames(object)
##	numberInformative <- list()
	##For each trio, fit the HMM
	for(i in seq(along=unique(familyId))){
		if(verbose) cat("Family ", unique(familyId)[i], ", ")
		fId <- unique(familyId)[i]
		trio <- which(fId == as.character(object$familyId))
		trioSet <- object[, trio]
		## TODO:
		## Remove the noinformative snps here.
		sampleNames(trioSet) <- trioSet$Sample.Name
		isBPI <- isBiparental.SnpCallSetPlus(trioSet, allowHetParent=hmmOptions[["allowHetParent"]])
		isInformative <- !is.na(isBPI)
		if(all(!isInformative)){
			fit[, i] <- 1
			next()
		}
		trioSet <- trioSet[isInformative, ]
		index <- match(featureNames(trioSet), rownames(fit))
		tau <- transitionProbability(chromosome=chromosome(trioSet),
					     position=position(trioSet),
					     TAUP=TAUP)
		emission <- computeBpiEmission.SnpCallSetPlus(trioSet, hmmOptions, isBPI=isBPI[isInformative])
		if(is.null(emission)) stop("not a father, mother, offspring trio")
		##NA's in the emission probabilities for SNPs that are not informative
		##noNAs <- rowSums(is.na(emission)) == 0
		##Just set the emission probability of noninformative SNPs to zero
		##emission[is.na(emission)] <- 0
##		if(any(noNAs)){ ## some of the markers are informative
		log.e <- array(log(emission), dim=c(nrow(trioSet), 1, 2), dimnames=list(featureNames(trioSet), fId, states))
		tmp <- as.matrix(as.integer(viterbi(initialStateProbs=log(initialP),
						    emission=log.e,
						    tau=tau[, "transitionPr"],
						    arm=tau[, "arm"],
						    normalIndex=1,
						    verbose=verbose)))
		##nInformative <- sum(isInformative & tmp != 1, na.rm=TRUE)
		## Chromosome |----------------------------------------------------------|
		##            |******|111**********111*1*1|******************************|
		## ->         |1111111111111111111111111111111111111111111111111111111111|
		if(all(tmp == 1, na.rm=TRUE)){
			fit[, i] <- 1
		} else {
			##            |******|111*********111*1*11|*******|2222****22|***********|
			## ->         |11111111111111111111111111111111111|2222222222|11111111111|
			##tmp does not have the asterisks
			##Calculate breaks on the tmp vector
			##For the altered breaks only, find the indices of all markers in the larger matrix.
			##Change the indices of all markers within a break to the altered state
			##Assign the remaining markers to the normal state.
			tmp2 <- breaks(tmp, states=states,
				       position=tau[, "position"],
				       chromosome=tau[, "chromosome"])
			tmp2 <- tmp2[tmp2$state=="notBPI", ]
			for(k in 1:nrow(tmp2)){
				xx <- tmp2[k, ]
				index <- which(position(object) >= xx[["start"]] & position(object) <= xx[["end"]])
				fit[index, i] <- 2
			}
			fit[is.na(fit[, i]), i] <- 1

		}
		##            |******|222**********222*2*2|******************************|
		## ->         |111111|22222222222222222222|111111111111111111111111111111|
	}
##	if(!all(is.na(fit))){
	tmp <- breaks(x=fit[, , drop=FALSE],
		      states=states,
		      position=position(object),
		      chromosome=chromosome(object))
	brks <- tmp[tmp$state == "notBPI" & !is.na(tmp$state), ]
	return(brks)
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



hmm.old <- function(states=c("BPI", "notBPI"),
			  initialP=c(0.99, 0.01),
			  TAUP=1e7,
			  prGtError=c(0.001, 0.01),
			  fns, verbose=FALSE,
			  samplesheet,
			  outdir){
	##Container for output
	brksAll <- vector("list", length(fns))

	##For each plate
	for(batch in seq(along=fns)){
		##find plate directory
		message("Batch ", batch, " of ", length(fns))
		PLATE <- names(fns)[batch]
		platedir <- file.path(outdir, PLATE)

		## Does plate directory exist.  If not, go to the next
		## plate.
		if(!file.exists(platedir)){
			message(platedir, " does not exist.")
			next()
		} else{
			##If there are fewer than 20 samples on the
			##plate, go to the next plate
			## *Note, for purposes of genotyping, plate is less important
			##  We should get crlmm genotypes for all the data
			if(length(list.files(platedir)) < 20){
				message(platedir, " not preprocessed.")
				next()
			}
		}

		## A list of breaks by chromosome
		brks <- vector("list", 22)
		for(CHR in 1:22){
			message("Chromosome ", CHR)

			## Load the chromosome output from crlmm
			fn <- paste("object_", CHR, ".rda", sep="")
			if(!file.exists(file.path(platedir, fn))){  ##check whether the file exists
				warning("file does not exist")
				next()
			} else{
				## We had problems loading a few of the samples
				message("Loading ", fn, "...")
				err <- tryCatch(load(file.path(platedir, fn)), error=function(e) return(NULL))
				if(is.null(err)){
					next("something wrong with the connection.  Need to preprocess/ genotype again.")
				}
			}
			## Add feature and sample annotation to the object object
			## * Would be better do do this during the estimation than post-hoc
			if(!("familyId" %in% varLabels(crlmmSetList[[1]]))){
				phenoData(crlmmSetList[[1]]) <- addSampleSheetAnnotation(crlmmSetList, samplesheet)
				pData(crlmmSetList[[1]])[, "familyId"] <- as.character(pData(crlmmSetList[[1]])[, "familyId"])
				pData(crlmmSetList[[1]])[, "familyMember"] <- as.character(pData(crlmmSetList[[1]])[, "familyMember"])
				featureData(crlmmSetList[[1]]) <- addFeatureAnnotation(crlmmSetList)
				crlmmSetList <- crlmmSetList[order(chromosome(crlmmSetList), position(crlmmSetList)), ]
				message("Saving crlmmSetList_", CHR, ".rda")
				save(crlmmSetList, file=file.path(platedir, fn))
			}
			tau <- transitionProbability(chromosome=chromosome(crlmmSetList),
						     position=position(crlmmSetList),
						     TAUP=TAUP)
			fit <- matrix(NA, nrow(crlmmSetList), ncol=length(unique(crlmmSetList[[1]]$familyId)))
			colnames(fit) <- unique(crlmmSetList[[1]]$familyId)

			##For each trio, fit the HMM
			for(i in seq(along=unique(crlmmSetList[[1]]$familyId))){
				if(ncol(crlmmSetList) < 3) next()
				id <- unique(crlmmSetList[[1]]$familyId)[i]
				trio <- which(id == as.character(crlmmSetList[[1]]$familyId))
				if(length(trio) != 3){
					if(verbose) message("Not three members in trio ", i, ".  Skipping...")
					next()
				}
				trioSet <- crlmmSetList[, trio]
				emission <- computeEmission(trioSet, prGtError)
				if(all(emission == 0)){
					fit[, i] <- rep(1, nrow(fit))
					next()  ## not a F, M, O trio
				}
				isInform <- rowSums(is.na(emission)) == 0
				emission[is.na(emission)] <- 0
				if(any(isInform)){
					e <- array(emission, dim=c(nrow(trioSet), 1, 2), dimnames=list(featureNames(trioSet), id, states))
					fit[, i] <- as.integer(viterbi(initialStateProbs=log(initialP),
								       emission=e,
								       tau=tau[, "transitionPr"],
								       arm=tau[, "arm"],
								       normalIndex=1, verbose=verbose))
				} else {
					message("Non-mendelian inheritance was not observed in any of the informative SNPs for family", id)
					fit[, i] <- rep(1, nrow(fit))
				}
			}
			if(!all(is.na(fit))){
				tmp <- breaks(x=fit[, , drop=FALSE],
					      states=states,
					      position=tau[, "position"],
					      chromosome=tau[, "chromosome"])
				brks[[CHR]] <- tmp[tmp$state == "notBPI" & !is.na(tmp$state), ]
			}
		}
		nullelements <- sapply(brks, is.null)
		brks <- brks[!nullelements]
		noHits <- sapply(brks, nrow) == 0
		brks <- brks[!noHits]
		if(length(brks) > 1){
			brks <- do.call(rbind, brks)
		} else brks <- brks[[1]]
		save(brks, file=file.path(platedir, "brks.rda"))
		brksAll[[batch]] <- brks
	}
}


getFamilyId <- function(object) substr(object$Sample.Name, 1, 5)
getIndividualId <- function(object) substr(object$Sample.Name, 7, 8)
getIndividualId2 <- function(object) substr(object$Sample.Name, 1, 8)


computeBpiEmission.SnpCallSetPlus <- function(object, hmmOptions, isBPI){
	states <- hmmOptions[["states"]]
	prGtError <- hmmOptions[["prGtError"]]
	emission <- matrix(NA, nrow(object), ncol=2)
	colnames(emission) <- states
	emission[isBPI==TRUE,  "BPI"] <-  1-prGtError["BPI"]
	emission[isBPI==FALSE, "BPI"] <- prGtError["BPI"] ##Mendelian inconsistancy
	emission[isBPI==TRUE,  "notBPI"] <- prGtError["notBPI"]   ## biparental inheritance, but true state is not Biparental
	emission[isBPI==FALSE, "notBPI"] <- 1-prGtError["notBPI"] ## Mendelian inconsistancy

##	##M/F AA, F/M BB, O AB
##	##isHet <- offspringHeterozygous(object)  ##offspring is heterozygous
##	biparental.0 <- F==1 & M == 3 & O == 2  #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	biparental.1 <- F==3 & M == 1 & O == 2  #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	emission[biparental.0 | biparental.1, "BPI"] <- 1-prGtError[1]
##	emission[biparental.0 | biparental.1, "notBPI"] <- prGtError[2]
##	##M/F AA, F/M BB, O AA or BB
##	biparental.2 <- F==1 & M == 3 & (O == 1 | O == 3) #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	biparental.3 <- F==3 & M == 1 & (O == 1 | O == 3) #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	emission[biparental.2 | biparental.3, "BPI"] <- prGtError[1]
##	emission[biparental.2 | biparental.3, "notBPI"] <- 1-prGtError[2]
##	##
##	## M/F AA, F/M BB, O AB
##	biparental.4 <- F == 1 & M == 2 & O == 2 #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	biparental.5 <- F == 2 & M == 1 & O == 2 #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	emission[biparental.4 | biparental.5, "BPI"] <- 1-prGtError[1]
##	emission[biparental.4 | biparental.5, "notBPI"] <- prGtError[2]
##	##
##	## F AA, M AB, O AB
##	## F AA, M AB, O BB
##	biparental.6 <- F == 2 & M == 1 & O == 3  #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	biparental.7 <- F == 1 & M == 2 & O == 3  #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	emission[biparental.6 | biparental.7, "BPI"] <- prGtError[1]
##	emission[biparental.6 | biparental.7, "notBPI"] <- 1-prGtError[2]
##	## M BB, F AB, O AB
##	biparental.10 <- F == 2 & M == 3 & O == 2 #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	## F BB, M AB, O AB
##	biparental.12 <- F == 3 & M == 2 & O == 2 #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
##	emission[biparental.10 | biparental.12, "BPI"] <- 1-prGtError[1]
##	emission[biparental.10 | biparental.12, "notBPI"] <- prGtError[2]
##	## M AA, F AB, O BB
##	## M BB, F AB, O AA
##	biparental.11 <- F == 2 & M == 3 & O == 1 #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	biparental.12 <- F == 3 & M == 2 & O == 1 #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
##	emission[biparental.11 | biparental.12, "BPI"] <- prGtError[1]
##	emission[biparental.11 | biparental.12, "notBPI"] <- 1-prGtError[2]
	return(emission)
}


computeEmission <- function(object, prGtError){
	sns <- sapply(object[[1]]$familyMember, who)
	whoisit <- sapply(object[[1]]$familyMember, who)
	father <- which(whoisit == "father")
	mother <- which(whoisit == "mother")
	offspring <- which(whoisit == "offspring")
	F <- calls(object[, father])
	M <- calls(object[, mother])
	O <- calls(object[, offspring])
	ncols <- c(ncol(F), ncol(M), ncol(O))
	if(any(ncols != 1)){
		message("Family ", unique(object[[1]]$familyId), " is not a Father, Mother, Offspring Trio.")
		emission <- matrix(0, nrow(object), ncol=2)
		colnames(emission) <- states
		return(emission)
	}
	emission <- matrix(NA, nrow(object), ncol=2)
	colnames(emission) <- states
	##M/F AA, F/M BB, O AB
	##isHet <- offspringHeterozygous(object)  ##offspring is heterozygous
	biparental.0 <- F==1 & M == 3 & O == 2  #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	biparental.1 <- F==3 & M == 1 & O == 2  #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	emission[biparental.0 | biparental.1, "BPI"] <- 1-prGtError[1]
	emission[biparental.0 | biparental.1, "notBPI"] <- prGtError[2]
	##M/F AA, F/M BB, O AA or BB
	biparental.2 <- F==1 & M == 3 & (O == 1 | O == 3) #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	biparental.3 <- F==3 & M == 1 & (O == 1 | O == 3) #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	emission[biparental.2 | biparental.3, "BPI"] <- prGtError[1]
	emission[biparental.2 | biparental.3, "notBPI"] <- 1-prGtError[2]
	##
	## M/F AA, F/M BB, O AB
	biparental.4 <- F == 1 & M == 2 & O == 2 #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	biparental.5 <- F == 2 & M == 1 & O == 2 #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	emission[biparental.4 | biparental.5, "BPI"] <- 1-prGtError[1]
	emission[biparental.4 | biparental.5, "notBPI"] <- prGtError[2]
	##
	## F AA, M AB, O AB
	## F AA, M AB, O BB
	biparental.6 <- F == 2 & M == 1 & O == 3  #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	biparental.7 <- F == 1 & M == 2 & O == 3  #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	emission[biparental.6 | biparental.7, "BPI"] <- prGtError[1]
	emission[biparental.6 | biparental.7, "notBPI"] <- 1-prGtError[2]
	## M BB, F AB, O AB
	biparental.10 <- F == 2 & M == 3 & O == 2 #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	## F BB, M AB, O AB
	biparental.12 <- F == 3 & M == 2 & O == 2 #Pr(O | biparental)=1-0.001, Pr(O | not biparental) = 0.01
	emission[biparental.10 | biparental.12, "BPI"] <- 1-prGtError[1]
	emission[biparental.10 | biparental.12, "notBPI"] <- prGtError[2]
	## M AA, F AB, O BB
	## M BB, F AB, O AA
	biparental.11 <- F == 2 & M == 3 & O == 1 #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	biparental.12 <- F == 3 & M == 2 & O == 1 #Pr(O | biparental)=0.001, Pr(O | not biparental) = 1-0.01
	emission[biparental.11 | biparental.12, "BPI"] <- prGtError[1]
	emission[biparental.11 | biparental.12, "notBPI"] <- 1-prGtError[2]
	return(emission)
}

isOverlap <- function(x, y){##x = bpi results, y= penn.cnv
	x.starts <- x[, "start"]
	x.ends <- x[, "end"]
	y.start <- y[["start"]]
	y.end <- y[["end"]]
	isXintervalInYinterval(x.starts, x.ends, y.start, y.end)
}

doesAnyXoverlapY <- function(x, y){
	if(nrow(y) > 1) stop("y should have one row (one segment)")
	isOverlap(x, y)
}

##Determine if x is in the interval y
whichPennCnvRows <- function(x, y){
	if(nrow(y) > 1) stop("should be 1 row in y data.frame")
	y.start <- y[["start"]]
	y.end <- y[["end"]]
	x.starts <- x[, "StartPosition"]
	x.ends <- x[, "EndPosition"]
	isXintervalInYinterval(x.starts, x.ends, y.start, y.end)
}

rerunBpiHMM <- function(notBPI){
	data(samplesheet, package="Beaty")
	samplesheet <- samplesheet[-grep("CIDR", samplesheet[, "Sample.Name"]), ]
	individualIds <- substr(samplesheet[, "Sample.Name"], 1, 8)
	plates <- unique(samplesheet[match(notBPI[, "ID"], individualIds), "Sample.Plate"])
	platedir <- file.path("/thumper/ctsa/beaty/scharpf/crlmmOut", plates)
	todo <- rep(NA, length(plates))
	for(i in seq(along=plates)){
		todo[i] <- ifelse(file.exists(file.path(platedir[i], "brks.rda")), FALSE, TRUE)
	}
	platedir[todo]
}

extractSNRandDates <- function(dirs, UNLINK){
	scanDateList <- snr <- list()
	for(i in seq(along=dirs)){
		cat(".")
		load(file.path(dirs[i], "crlmmSetList_22.rda"))
		callSetPlus <- as(crlmmSetList, "SnpCallSetPlus")
		if(!("familyId" %in% varLabels(callSetPlus))){
			phenoData(callSetPlus) <- addSampleSheet(callSetPlus)
		}
		snr[[i]] <- callSetPlus$SNR
		names(snr[[i]]) <- callSetPlus$Sample.Name
		##snr[[i]] <- crlmmSetList[[1]]@phenoData@data$SNR
		if(UNLINK){
			if(file.exists(file.path(dirs[i], "normalizedIntensities.rda")))
				unlink(file.path(dirs[i], "normalizedIntensities.rda"))
			if(file.exists(file.path(dirs[i], "snpsetObject.rda")))
				unlink(file.path(dirs[i], "snpsetObject.rda"))
		}
		##the protocol data should be added to the crlmmSetList object
		if(file.exists(file.path(dirs[i], "rgFile.rda"))){
			load(file.path(dirs[i], "rgFile.rda"))
			scanDates <- pData(protocolData(RG))
			scanDateList[[i]] <- scanDates[match(rownames(pData(crlmmSetList[[1]]@phenoData)), rownames(scanDates)), ]
		} else{
			scanDateList[[i]] <- NULL
		}
	}
	names(snr) <- dirs
	##snr <- unlist(snr)
	message("Saving SNR and scanDateList to ", dirs[1])
	save(snr, file=file.path(dirname(dirs[1]), "snr.rda"))
	save(scanDateList, file=file.path(dirname(dirs[1]), "scanDateList.rda"))
}

getPennSegments <- function(datadir, deletions.only=TRUE){
	if(file.exists("~/projects/Beaty/data/penn.seg.rda")){
		message("Loading penn.seg.rda")
		load("~/projects/Beaty/data/penn.seg.rda")
		penn.seg <- get("penn.seg")
		if(deletions.only)
			penn.seg <- penn.seg[sapply(penn.seg[, "TrioState"], isDeletion), ]
		return(penn.seg)
	}
	fns <- list.files(datadir, full.names=TRUE)
	penn.seg <- list()
	##fn <- substr(basename(fns), 1, 5)
	for(i in seq(along=fns)){
		cat(".")
		penn.seg[[i]] <- read.delim(fns[i], as.is=TRUE)
		##penn.seg[[i]]$fn <- as.character(fn[i])
		if(i == length(fns)) cat("\n")
	}
	penn.seg <- do.call("rbind", penn.seg)
	## The ID column in penn.seg looks like XXXXX_YY
	## YY = 01: offspring
	## YY = 02: mother
	## YY = 03: father
	##---------------------------------------------------------------------------
	## Exclusions:
	## Drop names that begin with CIDR
	data(samplesheet, package="Beaty")
	samplesheet <- samplesheet[-grep("CIDR", samplesheet[, "Sample.Name"]), ]
	## Remove the offspring that were not part of a trio.
	familyId <- substr(samplesheet[, "Sample.Name"], 1, 5)
	individualId <- substr(samplesheet[, "Sample.Name"], 1, 8)
	individualId.family <- split(individualId, familyId)
	isCompleteTrio <- function(x){
		length(x) >= 3 & length(grep("_01", x)) == 1 & length(grep("_02", x)) == 1 & length(grep("_03", x))
	}
	completeTrios <- sapply(individualId.family, isCompleteTrio)
	individualId.family <- individualId.family[completeTrios]
	individualId <- unlist(individualId.family)
					# remove 60 individuals that did not have any alterations identified by penncnv
	individualId <- individualId[individualId %in% penn.seg[, "ID"]]
	##Remove individuals from penn cnv output that are not a member of a complete trio
	penn.seg <- penn.seg[penn.seg[, "ID"] %in% individualId, ]
	## Finally, remove the non-offspring samples  (do this last)
	isOffspring <- substr(penn.seg$ID, 7,8) == "01"
	penn.seg <- penn.seg[isOffspring, ]
	colnames(penn.seg)[c(1:4, 6)] <- c("chr", "start", "end", "nprobes", "sample")
	save(penn.seg, file="~/projects/Beaty/data/penn.seg.rda")

	if(deletions.only) penn.seg <- penn.seg[sapply(penn.seg[, "TrioState"], isDeletion), ]
	return(penn.seg)
}

getDirs <- function(path="/thumper/ctsa/beaty/scharpf/crlmmOut",
		    CIDR=FALSE,
		    MIN.SAMPLES=3){
	data(samplesheet, package="Beaty")
	if(!CIDR){
		message("Excluding CIDR samples")
		samplesheet <- samplesheet[-grep("CIDR", samplesheet[, "Sample.Name"]), ]
	}
	fns <- split(samplesheet$filenames, samplesheet[, "Sample.Plate"])
	dirs <- names(fns)[sapply(fns, length) > MIN.SAMPLES]
	file.path(path, dirs)
}

getBpiHmm <- function(filename="brks.rda", ...){
	dirs <- getDirs(...)
	brksAll <- list(); k <- 1
	for(i in seq(along=dirs)){
		if(file.exists(file.path(dirs[i], filename))){
			load(file.path(dirs[i], filename))
		} else next()
		brksAll[[k]] <- brks
		k <- k+1
	}
	brksAll <- do.call("rbind", brksAll)
	brks <- brksAll; rm(brksAll); gc()
	## The sample column in bpi.seg looks like XXXXX.
	## Add the _01 extension so that it corresponds to pennCnv
	brks[, "sample"] <- paste(brks[, "sample"], "_01", sep="")
	return(brks)
}

pennRegionsInBPI <- function(penn.seg, bpi.seg, load.it=TRUE){
	if(load.it & file.exists("~/projects/Beaty/inst/scripts/inBPI.rda")){
		message("Loading inBPI.rda")
		load("~/projects/Beaty/inst/scripts/inBPI.rda")
		return(inBPI)
	}
	inBPI <- rep(NA, nrow(penn.seg))
	penn.seg <- penn.seg[order(penn.seg$NumberSNPs, decreasing=TRUE), ]
	for(i in 1:nrow(penn.seg)){
		if(i %% 100 == 0) cat(".")
		##for each alteration identified by penn cnv, assess whether the region is covered by the BPI HMM
		y <- penn.seg[i, ]
		j <- grep(y[["ID"]], bpi.seg[, "sample"])
		if(length(j) < 1){
			inBPI[i] <- FALSE
			next()
		} else{## determine if there is any overlap
			##Check if on the same chromosome
			x <- bpi.seg[j, , drop=FALSE]
			if(y[["Chromosome"]] %in% x[, "chr"]){
				tmp <- isOverlap(x=bpi.seg[j, , drop=FALSE], y=y)
				inBPI[i] <- any(tmp) ##any overlap
			} else {
				##regions not on same chromosome
				inBPI[i] <- FALSE
			}
		}
	}
	save(inBPI, file="inBPI.rda")
	return(inBPI)
}

match.data.frames <- function(x, table){
	indicator <- rep(NA, nrow(x))
	for(i in 1:nrow(x)){
		if(i %% 100 == 0) cat(".")
		##for each alteration identified in x, assess whether the region is covered in 'table'
		xi <- x[i, ]
		##			  y <- penn.seg[i, ]
		##			  j <- grep(y[["ID"]], bpi.seg[, "sample"])
		j <- grep(xi[["sample"]], table[, "sample"])
		if(length(j) < 1){
			indicator[i] <- FALSE
			next()
		} else{## determine if there is any overlap
			##Check if on the same chromosome
			tablej <- table[j, , drop=FALSE]
			##x <- bpi.seg[j, , drop=FALSE]
			if(xi[["chr"]] %in% tablej[, "chr"]){
				##				  if(y[["Chromosome"]] %in% x[, "chr"]){
				tmp <- doesAnyXoverlapY(x=tablej, y=xi)
				##tmp <- isOverlap(x=bpi.seg[j, , drop=FALSE], y=y)
				indicator[i] <- any(tmp) ##any overlap
			} else {
				##regions not on same chromosome
				indicator[i] <- FALSE
			}
		}
	}
	indicator
}

setMethod("%in%", signature(x="data.frame", table="data.frame"),
	  function(x, table){
		  match.data.frames(x, table)
	  })


bpiSubmitter <- function(dirs.split, sleep=0, allowHetParent, version=1, splitSize=5){
	if(missing(allowHetParent)) stop("must specify allowHetParent")
	for(j in seq(along=dirs.split)){
		sink("temp")
		cat("j <- ", j, "\n")
		cat("allowHetParent <- ", allowHetParent, "\n")
		cat("version <- ", version, "\n")
		cat("splitSize <- ", splitSize, "\n")
		sink()
		fn <- paste("bpi_", j, ".R", sep="")
		if(file.exists(fn)) system(paste("rm", fn))
		system(paste("cat temp bpiScript.R >", fn))
		##Requires very little RAM (700MB)
		system(paste("cluster 2G", fn))
		Sys.sleep(sleep)
	}
}

crlmmCalls_submit <- function(){
	for(k in 2:11){
		sink("temp")
		cat("k <- ", k, "\n")
		sink()
		fn <- paste("crlmmGt_", k, ".R", sep="")
		if(file.exists(fn)) system(paste("rm", fn))
		system(paste("cat temp crlmmGenotypes.R >", fn))
		##Requires very little RAM (700MB)
		system(paste("cluster 4G", fn))
	}
}


hmm.wrapper <- function(dirs, chromosomes=1:22, hmmOptions, filename="brks"){
	brksAll <- vector("list", length(dirs))
	for(i in seq(along=dirs)){
		brks <- vector("list", length(chromosomes))
		write.header <- TRUE
		k <- 1
		for(CHR in chromosomes){
			cat("\nChromosome ", CHR, "\n")
			##this is just one chromosome
			load(file.path(dirs[i], paste("crlmmSetList_", CHR, ".rda", sep="")))
			crlmmSetList <- get("crlmmSetList")
			callSetPlus <- as(crlmmSetList, "SnpCallSetPlus")
			rm(crlmmSetList); gc()
			if(!("chromosome" %in% fvarLabels(callSetPlus))){
				featureData(callSetPlus) <- addFeatureAnnotation(callSetPlus)
			}
			if(!("familyId" %in% varLabels(callSetPlus))){
				phenoData(callSetPlus) <- addSampleSheet(callSetPlus)
			}
			if(!("isSnp" %in% fvarLabels(callSetPlus))){
				fData(callSetPlus)$isSnp <- isSnp(callSetPlus)
			}
			callSetPlus <- callSetPlus[, -grep("CIDR", callSetPlus$Sample.Name)]
			pData(callSetPlus)$familyMember <- getIndividualId(callSetPlus)
			callSetPlus <- callSetPlus[isSnp(callSetPlus), ]
			callSetPlus <- callSetPlus[order(position(callSetPlus)), ]
			brks[[k]] <- hmm.SnpCallSetPlus(object=callSetPlus, hmmOptions=hmmOptions)
			if(is.null(brks[[k]])) break() ## Go to the next plate
			rm(callSetPlus); gc()
			k <- k+1
		}
		backup <- brks
		brks <- do.call("rbind", brks)
		save(brks, file=file.path(dirs[i], paste(filename, ".rda", sep="")))
	}
	return(brks)
}

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


plotCandidates <- function(){


}


getTdtSnps <- function(x){
	chrom <- unique(x[, "CHR"])
	dirs <- getDirs(CIDR=FALSE, MIN.SAMPLES=3)
	aa <- list()
	bb <- list()
	gg <- list()

	for(i in seq(along=chrom)){
		CHR <- chrom[i]
		xx <- x[x[, "CHR"] == CHR, ]
		id <- xx[, "SNP"]
		a <- list()
		b <- list()
		gt <- list()
		for(p in seq(along=dirs)){
			message("Loading:", file.path(dirs[p], paste("crlmmSetList_", CHR, ".rda", sep="")))
			load(file.path(dirs[p], paste("crlmmSetList_", CHR, ".rda", sep="")))
			callSetPlus <- as(crlmmSetList, "SnpCallSetPlus")
			if(!"familyId" %in% varLabels(callSetPlus)){
				phenoData(callSetPlus) <- addSampleSheet(callSetPlus)
			}
			callSetPlus <- callSetPlus[, -grep("CIDR", callSetPlus$Sample.Name)]
			a[[p]] <- A(callSetPlus)[match(id, featureNames(callSetPlus)), ]
			b[[p]] <- B(callSetPlus)[match(id, featureNames(callSetPlus)), ]
			gt[[p]] <- calls(callSetPlus)[match(id, featureNames(callSetPlus)), ]
		}
		aa[[i]] <- do.call("cbind", a)
		bb[[i]] <- do.call("cbind", b)
		gg[[i]] <- do.call("cbind", gt)
		rm(a, b, gt)
	}
	results <- list(aa=aa, bb=bb, gg=gg)
	save(results, file="results.rda")
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

isHomozygousDeletionInOffspring <- function(x){
	if(length(grep("-", x)) > 0){
		tmp <- strsplit(x, "_")[[1]]
		state <- substr(tmp, 3, 3)
		state <- ifelse(any(state == 1), TRUE, FALSE)
	} else{
		state <- as.integer(substr(x, 3, 3))
		state <- ifelse(state == 1, TRUE, FALSE)
	}
	state
}

##getInfMarkersPennSeg <- function(segments, hmmOptions, index=1:500, version=1, save.it=TRUE){
numberInfMarkers <- function(segments, hmmOptions, index=1:500, version=1, save.it=TRUE){
	if(missing(segments)) stop("must provide PennCNV results")
	if(missing(hmmOptions)) stop("must provide hmm options")
	BPI <- matrix(NA, nrow(segments), 3)
	rownames(BPI) <- rownames(segments)
	colnames(BPI) <- c("nBiparental", "nNotBiparental", "nNonInformative")
	data(samplesheet)
	samplesheet <- samplesheet[-grep("CIDR", samplesheet[, "Sample.Name"]), ]
	path <- "/thumper/ctsa/beaty/scharpf/crlmmOut"
	x <- segments[1, ]
	id <- x[["sample"]]
	chr <- x[["chr"]]
	plate <- samplesheet[match(id, substr(samplesheet[, "Sample.Name"], 1, 8)), "Sample.Plate"]
	platedir <- file.path(path, plate)
	fn <- file.path(platedir, paste("crlmmSetList_", chr, ".rda", sep=""))
	message("Loading ", basename(fn))
	load(fn)
	crlmmSetList <- get("crlmmSetList")
	##Convert this to the current class definition
	callSetPlus <- as(crlmmSetList, "SnpCallSetPlus")
	## Extract the trio in questions
	if(!("familyId" %in% varLabels(callSetPlus))){
		phenoData(callSetPlus) <- addSampleSheet(callSetPlus)
	}
	callSetPlus <- callSetPlus[, -grep("CIDR", callSetPlus$Sample.Name)]
	if(!("chromosome" %in% fvarLabels(callSetPlus))){
		featureData(callSetPlus) <- addFeatureAnnotation(callSetPlus)
	}
	callSetPlus$familyId <- substr(callSetPlus$Sample.Name, 1, 5)
	pData(callSetPlus)$familyMember <- getIndividualId(callSetPlus)
	familyId <- unique(substr(segments[, "sample"], 1, 5))
	for(i in seq(along=familyId)){
		if(hmmOptions[["verbose"]]) cat(familyId[i], "\n")
		trioSet <- callSetPlus[, grep(familyId[i], callSetPlus$familyId)]
		if(ncol(trioSet) < 3){
			J <- grep(familyId[i], substr(segments[, "sample"], 1, 5))
			BPI[J, ] <- NA
			next()  ## go to next family
		}
		isBpi <- isBiparental.SnpCallSetPlus(trioSet, allowHetParent=hmmOptions[["allowHetParent"]])
		##how many rows are for this family
		J <- grep(familyId[i], substr(segments[, "sample"], 1, 5))
		for(j in J){
			region <- position(trioSet) >= segments[j, "start"] & position(trioSet) <= segments[j, "end"]
			if(length(j) < 1){
				BPI[j, ] <- NA
				next() ##go to next segment
			}
			nNonInformative <- sum(is.na(isBpi[region]))
			nInformative <- sum(!is.na(isBpi[region]))
			nNotBpi <- sum(isBpi[region] == FALSE, na.rm=TRUE)
			nBpi <- sum(isBpi[region] == TRUE, na.rm=TRUE)
			stopifnot((nNotBpi+nBpi) == nInformative)
			BPI[j, "nBiparental"] <- nBpi
			BPI[j, "nNotBiparental"] <- nNotBpi
			BPI[j, "nNonInformative"] <- nNonInformative
		}
	}
	return(BPI)
}

calculateNumberInformative <- function(outfile.csv, job, penn.seg){
	if(missing(job)) stop("must specify job")
	message("This will fail if you are not on a fat node!")
	message("Loading the penn segmentation for the deletions only")
	if(missing(penn.seg)){
		penn.seg <- getPennSegments(datadir="/thumper/ctsa/beaty/holger/penncnv/jointDat",
					    deletions.only=TRUE)
		message("Creating rowlabels for the penn segmentation that will be used to label the rows in the output.csv")
		rownames(penn.seg) <- paste("segment", 1:nrow(penn.seg), sep="_")
	}
	hmmOpts <- hmm.options(allowHetParent=FALSE)
	hmmOpts[["verbose"]] <- FALSE

	## Processed data is stored by plate and chromosome.  Create a factor
	## so that the data for all segments on a given plate and chromosome
	## is only read once.
	data(samplesheet)
	samplesheet <- samplesheet[-grep("CIDR", samplesheet[, "Sample.Name"]), ]
	##Split rows by chromosome and plate.  Process all of the rows
	plateChrom.factor <- paste(samplesheet[match(penn.seg[, "sample"], substr(samplesheet[, "Sample.Name"], 1, 8)), "Sample.Plate"], penn.seg[, "chr"], sep="_")
	index.plateChrom <- split(1:nrow(penn.seg), plateChrom.factor)
	L <- length(index.plateChrom)
	##split into 10 jobs
	index.jobs <- split(seq(along=index.plateChrom), rep(1:10, each=L/10, length.out=L))
	index <- index.plateChrom[index.jobs[[job]]]
	for(j in seq(along=index)){
		##for(i in 1:441){
		cat(j, "\n")
		rows <- index[[j]]
		tmp <- penn.seg[rows, ]
		X <- numberInfMarkers(tmp, hmmOpts, 1:nrow(tmp), version=1)
		write.csv(X, file=outfile.csv, append=TRUE, quote=FALSE, col.names=TRUE)
	}
}


tdtMatrix <- function(){
	tmp <- matrix(NA, nrow=3^3, ncol=7)
	colnames(tmp) <- c("father", "mother", "offspring", "transmittedB", "untransmittedB", "transmittedA", "untransmittedA")
	counter <- 1
	for(i in c("AA", "AB", "BB")){##genotype of father
		for(j in c("AA", "AB", "BB")){ ## genotype of mother
			for(k in c("AA", "AB", "BB")){ ##genotype of offspring
				tmp[counter, 1:3] <- c(i, j, k)
				counter <- counter+1
			}
		}
	}
	tmp[, "transmittedB"] <-   c(0, 0, 0, 0, 1, 1, 0, 1, 2, 0, 1,1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2)
	tmp[, "untransmittedB"] <- c(0, 0, 0, 1, 0, 0, 2, 1, 0, 1, 0, 0, 2, 1, 0, 3, 2, 1, 2, 1, 0, 3, 2, 1, 4, 3, 2)
	tmp[, c(6, 7)] <-   matrix(c(2, 2,
				     1, 3,
				     0, 4,
				     2, 1,
				     1, 1,
				     0, 3,
				     2, 0,
				     1, 1,
				     0, 2,
				     2, 1,
				     1, 2,
				     0, 3,
				     2, 0,
				     1, 1,
				     0, 2,
				     1, 0,
				     1, 0,
				     0, 1,
				     2, 0,
				     1, 1,
				     0, 2,
				     1, 0,
				     1, 0,
				     0, 1,
				     0, 0,
				     0, 0,
				     0, 0), byrow=TRUE, ncol=2)
	return(tmp)
}


computeTdt <- function(crlmmCalls){
	##match(colnames(crlmmCalls), substr(samplesheet[, "Sample.Name"], 1, 8))
	col.index <- split(1:ncol(crlmmCalls), substr(colnames(crlmmCalls), 1, 5))
	tdtM <- tdtMatrix()
	TU <- vector("list", length(col.index))
	for(m in seq(along=col.index)){
		if(length(col.index[[m]]) < 3) next()  ##next family
		gts <- crlmmCalls[, col.index[[m]]]
		colnames(gts) <- who3(colnames(gts))
		gts <- gts[, match(c("father", "mother", "offspring"), colnames(gts))]
		##possibilities
		tmp <- matrix(NA, nrow(crlmmCalls), 2)
		colnames(tmp) <- c("T", "U")
		counter <- 1
		for(i in 1:3){##genotype of father
			for(j in 1:3){ ## genotype of mother
				for(k in 1:3){ ##genotype of offspring
					index <- which(gts[, 1] == i & gts[, 2] == j & gts[, 3] == k)
					tmp[index, "T"] <- as.integer(tdtM[counter, "transmittedB"]) + as.integer(tdtM[counter, "transmittedA"])
					tmp[index, "U"] <- as.integer(tdtM[counter, "untransmittedB"]) + as.integer(tdtM[counter, "untransmittedA"])
					counter <- counter+1
				}
			}
		}
		TU[[m]] <- tmp
	}
	matFun <- function(lis, FUN){
		lis <- lis[!sapply(lis, is.null)]
		if(!is.list(lis) || !all(sapply(lis, is.matrix)))
			stop("'lis' must be a list containing 2-dimensional arrays")
		dims <- sapply(lis, dim)
		n <- dims[1, 1]
		p <- dims[2, 1]
		if(!all(n == dims[1, ]) || !all(p == dims[2, ]))
			stop("the matrices must have the same dimensions")
		mat <- matrix(unlist(lis), n * p, length(lis))
		matrix(apply(mat, 1, FUN), n, p)
	}
	err <- tryCatch(tmp <- matFun(TU, "sum"), error=function(e) NULL)
	if(!is.null(err)){
		TU <- tmp
	} ##otherwise, keep TU  in list form
	return(TU)
}


my.qq.func=function(p,f,pval=TRUE,gc=FALSE,tn,hc=0.99,hm=100,...){
  if(missing(f)) f=rep(0,length(p))
  wh=which(is.na(p)|is.na(f))
  if(length(wh)>0){
    p=p[-wh]
    f=f[-wh]
  }
  wh=order(p)
  p=p[wh]
  f=f[wh]
  if(gc){
    z=qchisq(p,1,lower=F)
    gcp=median(z)/qchisq(0.5,1)
    cat("The genomic control parameter is ",gcp,"\n")
    z=z/gcp
    pchisq(z,1,lower=F)
  }
  n=length(p)
  x1=1:n
  x2=n+1-x1
  x=x1/(n+1)
  up=qbeta(0.975,x1,x2)
  lo=qbeta(0.025,x1,x2)
  if(pval==TRUE){
    uu=up
    up=-log10(lo)
    lo=-log10(uu)
    x=-log10(x)
    z=-log10(p)
  }
  else{
    up=qchisq(up,1)
    lo=qchisq(lo,1)
    x=qchisq(x,1)
    z=qchisq(p,1,lower=F)
  }
  z=rev(z)
  f=rev(f)
  x=rev(x)
  up=rev(up)
  lo=rev(lo)
  tt=NULL
  if(!missing(tn)){
    mx=10^tn
    tt=list(tn=mx,k=rev(x)[mx])
  }
  if(hm>1){
    nc=round(n*hc)
    wh=seq(1,nc,hm)
    wh=c(wh,(nc+1):n)
    z=z[wh]
    f=f[wh]
    x=x[wh]
    up=up[wh]
    lo=lo[wh]
  }
  return(list(z=z,f=f,x=x,up=up,lo=lo,tt=tt))
}

my.qq.plot=function(zz,cut,rmx=0,mt="",mt.cex=1,mt.line=NA,mgp=c(3,1,0),tn.cex=1,plab=T,xr,yr,...){
  if(!missing(cut)) zz$z[zz$z>cut]=cut
  if(missing(xr)) xr=c(0,1.02*max(zz$x))
  if(missing(yr)) yr=c(0,1.02*max(c(zz$up,zz$z,rmx)))
  if(plab){
    xl=expression(paste("expected  ",-log[10]," (p-value)",sep=""))
    yl=expression(paste("observed  ",-log[10]," (p-value)",sep=""))
  }
  else{
    xl=""
    yl=""
  }
  par(las=1)
  print(xr)
  plot(range(zz$x),range(c(0,zz$up)),type="n",xlim=xr,ylim=yr,xlab=xl,ylab=yl,...)
  axis(2,0:ceiling(yr)[2],...)
  polygon(c(zz$x,rev(zz$x)),c(zz$lo,rev(zz$up)),col="lightgrey",border=F)
  lines(c(0,max(zz$x)),c(0,max(zz$x)))
  cls=c("blue","red")
  cls=cls[zz$f+1]
  points(zz$x,zz$z,pch=20,col=cls,cex=0.5)
  par(mgp=mgp)
  if(length(zz$tt)==2){
    options(scipen=7)
    axis(3,zz$tt$k,zz$tt$tn,cex.axis=tn.cex)
    options(scipen=0)
  }
  title(mt,cex.main=mt.cex,line=mt.line)
  par(mgp=c(3,1,0))
}


loadCrlmmGenotypes <- function(){
	data(samplesheet)
	load("results.rda")
	a <- results[["aa"]]
	A <- do.call("rbind", a)
	b <- results[["bb"]]
	B <- do.call("rbind", b)
	g <- results[["gg"]]
	G <- do.call("rbind", g)
	sns <- substr(samplesheet[match(colnames(A), samplesheet[, "sampleNames"]), "Sample.Name"], 1, 8)
	colnames(A) <- colnames(B) <- colnames(G) <- sns

	individualId <- substr(colnames(G), 7, 8)
	individualId[individualId == "01"] <- "O" ##offspring
	individualId[individualId == "02"] <- "M" ##mother
	individualId[individualId == "02"] <- "F" ##father
	G <- G[, individualId != "04" & individualId != "05" & individualId != "06"]
	familyId <- substr(colnames(G), 1, 5)
	tmp <- split(1:ncol(G), familyId)
	##remove indices that are not trios
	tmp <- unlist(tmp[sapply(tmp, length) < 3])
	G <- G[, -tmp]

	A <- A[, match(colnames(G), colnames(A))]
	B <- B[, match(colnames(G), colnames(B))]

	pD <- annotatedDataFrameFrom(A, byrow=FALSE)
	fD <- annotatedDataFrameFrom(A, byrow=TRUE)

	callProbability <- matrix(NA, nrow(G), ncol(G))
	dimnames(callProbability) <- dimnames(G)
	callSet <- new("SnpCallSetPlus",
		       senseThetaA=A,
		       senseThetaB=B,
		       call=G,
		       callProbability=callProbability,
		       phenoData=pD,
		       featureData=fD,
		       annotation="human610quadv1b")
	return(callSet)
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

discordantCalls <- function(crlmmCallSet, BS, N){
	dns <- dimnames(BS)
	BS <- as.integer(BS)
	GG <- as.integer(calls(crlmmCallSet))
	notEqual <- rep(NA, length(BS))
	notEqual[!is.na(BS)] <- matrix(BS[!is.na(BS)] != calls(crlmmCallSet)[!is.na(BS)])
	notEqual <- matrix(notEqual, nrow(crlmmCallSet), ncol(crlmmCallSet))
	nE <- rowSums(notEqual, na.rm=TRUE)
	index <- which(nE > N)
	##BS <- matrix(BS, nrow(crlmmCallSet), ncol(crlmmCallSet))
	index
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

scaleSnr <- function(snr, min.scale, max.scale){
	tmp <- (snr-min(snr))/(max(snr)-min(snr))  ## 0 -> 1
	##max.scale <- 1
	##min.scale <- 0
	b <- 1/(max.scale-min.scale)
	a <- min.scale*b
	bg.scale <- (tmp + a)/b
	return(bg.scale)
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

addPhenoData <- function(fns){
	fns <- fns[-c(grep("23", fns), grep("24", fns))]
	##ped <- hapmapPedFile()  ##where is this function???
	ped <- read.csv("~/projects/Beaty/inst/extdata/HapMap_samples.csv")
	for(i in seq(along=fns)){
		cat(".")
		load(fns[i])
		fileExt <- strsplit(fns[i], "_")[[1]][2]
		sns <- sapply(sampleNames(cnSet), function(x) strsplit(x, "_")[[1]][[1]])
		ped <- ped[ped[, "coriellId"] %in% sns, ]
		stopifnot(nrow(ped) == ncol(cnSet))
		ped <- ped[match(sns, ped[, "coriellId"]), ]
		stopifnot(identical(ped[, "coriellId"], as.character(sns)))
		rownames(ped) <- sampleNames(cnSet)
		pD <- pData(cnSet)[, 1:3]
		pD2 <- cbind(pD, ped)
		pData(cnSet) <- pD2
		save(cnSet, file=fns[i])
	}
	if(file.exists(".RData")) unlink(".RData")

##	pD <- pData(cnSet)
##	pD2 <- pData(cnSet21)
##	pD <- cbind(pD2, pD[,4:8])
##	pD <- new("AnnotatedDataFrame", data=data.frame(pD), varMetadata=data.frame(labelDescription=colnames(pD), row.names=colnames(pD)))
##	phenoData(cnSet21) <- pD
##	cnSet <- cnSet21
##	save(cnSet, file=file.path(cnOpts[["outdir"]], "cnSet_21.rda"))
}

trio.fit <- function(cnSetFile){
	##i <- 10
	##chr <- as.numeric(sapply(basename(cnSetFiles), function(x) strsplit(strsplit(x, "_")[[1]][2], ".rda")[[1]][1])[i])
	chr <- strsplit(strsplit(basename(cnSetFile), "_")[[1]][2], ".rda")[[1]][1]
	if(chr > 22) next()
	if(chr == 21) next()
	cat("chromosome ", chr, "\n")
	load(cnSetFile)
	cnSet <- get("cnSet")
	cnSet$batch <- substr(sampleNames(cnSet), 13, 13)
	##trace(hmmOptions, browser)
	trioOpts <- hmmOptions(cnSet,
			       verbose=TRUE,
			       TAUP=1e10,
			       normalIndex=1,
			       altered2normal=0.5,
			       normal2altered=0.1,
			       trioHmm=TRUE)
	trios <- trioOpts[["trios"]]
	log.emission <- list()
	for(j in 1:nrow(trios)){
		cat(".")
		trioSet <- cnSet[, match(trios[j, ], sampleNames(cnSet))]
		isBPI <- isBiparental.SnpSuperSet(trioSet)
		isInformative <- !is.na(isBPI)
		if(all(!isInformative)){ log.emission[[j]] <- NULL; next()}
		trioSet <- trioSet[isInformative, ]
		isBPI <- isBPI[isInformative]
		log.emission[[j]] <- computeBpiEmission.SnpSuperSet(trioSet, trioOpts, isBPI=isBPI)
		rownames(log.emission[[j]]) <- featureNames(trioSet)
	}
	names(log.emission) <- trios[, "offspring"]
	rD <- list()
	for(j in seq(along=log.emission)){
		cat(".")
		trioOpts[["log.emission"]] <- array(log.emission[[j]], dim=c(nrow(log.emission[[j]]), 1, 2),
						    dimnames=list(rownames(log.emission[[j]]),
						    trios[j, "offspring"],
						    c("BPI", "notBPI")))
		object <- cnSet[match(rownames(log.emission[[j]]), featureNames(cnSet)), ]
		object <- object[, match(trios[j, ], sampleNames(object))]
		##trace(viterbi, browser)
		rD[[j]] <- hmm(object, trioOpts)
	}
	rD <- do.call(c, rD)
	## Do these regions appear plausible for de-novo events?
	rD <- rD[rD$LLR > 0, ]
	colnames(rD)[3] <- "informativeMarkers"
	rD <- rD[rD$informativeMarkers >= 3, ]
	## Find the total number of markers in each region
	## Exclude regions in which the percentage of informative markers is less than 2%
	rD$percentInformative <- rD$totalMarkers <- rep(NA, nrow(rD))
	for(k in 1:nrow(rD)){
		rD$totalMarkers[k] <- sum(position(cnSet) >= start(rD)[k] & position(cnSet) <= end(rD)[k])
	}
	rD$percentInformative <- rD$informativeMarkers/rD$totalMarkers
	rD <- rD[rD$percentInformative >= 0.02, ]
	return(rD)
}

combineElements <- function(rangedData){
	isnull <- sapply(rangedData, is.null)
	if(any(isnull))
		rangedData <- rangedData[!isnull]
	sampleId <- unlist(lapply(rangedData, function(x) x$sampleId))
	state <- unlist(lapply(rangedData, function(x) x$state))
	infM <- unlist(lapply(rangedData, function(x) x$informativeMarkers))
	totM <- unlist(lapply(rangedData, function(x) x$totalMarkers))
	perI <- unlist(lapply(rangedData, function(x) x$percentInformative))
	LLR <- unlist(lapply(rangedData, function(x) x$LLR))
	start <- unlist(lapply(rangedData, start))
	end <- unlist(lapply(rangedData, end))
	chr <- unlist(lapply(rangedData, space))
	rD <- RangedData(IRanges(start=start, end=end), space=chr,
			 sampleId=sampleId,
			 state=state,
			 informativeMarkers=infM,
			 totalMarkers=totM,
			 percentInformative=perI,
			 LLR=LLR)
	return(rD)
}

getHapmapTrios <- function(outdir){
	if(!file.exists("~/projects/Beaty/inst/extdata/trios_hapmap.rda")){
		load(file.path(outdir, "cnSet_20.rda"))
		object <- cnSet
		offspringId <- sampleNames(object)[object$fatherId != 0 & object$motherId != 0]
		trios <- as.matrix(t(sapply(offspringId, findFatherMother, object=object)))
		trios <- trios[rowSums(is.na(trios)) == 0, , drop=FALSE]
		colnames(trios) <- c("father", "mother", "offspring")
		trios
		save(trios, file="~/projects/Beaty/inst/extdata/trios_hapmap.rda")
	} else {
		load("~/projects/Beaty/inst/extdata/trios_hapmap.rda")
	}
	trios
}

getFamilyInfo <- function(phenoData){
	pheno <- read.csv2("~/projects/Beaty/inst/extdata/may_peds.csv", as.is=TRUE)
	merged <- merge(phenoData, pheno, by.x="CIDR_Name", by.y="cidr_name", all.x=TRUE)
	stopifnot(identical(phenoData$CIDR_Name, merged$CIDR_Name))
	rownames(merged) <- rownames(phenoData)
	return(merged)
}

constructPedigreeFromBsSet <- function(object){
	pedigree <- pData(object)[, c(5, 34:37)]
	pedigree$Abbrv.Name <- substr(pedigree$Sample.Name, 1, 8)
	save(pedigree, file="~/projects/Beaty/data/pedigree.rda")
}

getTriosIndex <- function(object){
	trios.index <- split(1:ncol(object), object$family)
	individ <- split(object$individ, object$family)
	trioid <- split(object$Sample.Name, object$family)
	trioid <- unlist(trioid)[unlist(individ) == 1]
	exclude <- sapply(individ, function(x) all(x != 1))
	trios.index <- trios.index[!exclude]
	trios.index <- trios.index[sapply(trios.index, length) >= 3]
	return(trios.index)
}

computeNumberDenovo <- function(object){
	x <- position(object)
	events <- matrix(NA, ncol(object), 3)
	for(j in 1:ncol(object)){
		if(j %% 100 == 0) cat(".")
		tmp <- D[, j]
		if(all(tmp == 0)) next()
		cs <- c(0, cumsum(diff(tmp) != 0))
		number.events <- split(tmp, cs)
		## the events are when the indicator is 1
		number.events <- length(number.events[sapply(number.events, unique) == 1])
		number.markers.per.event <- sapply(split(cs[-1], cs[-1]), length)
		median.number.markers.per.event <- median(number.markers.per.event)
		## a few events have just one marker.  Why?
		## cbs is run separately on each individual in the trio.  A length of 1 just means that the intersection had length 1.
		##indices.one.marker <- split(1:nrow(object), cs)
		##	indices.one.marker <- indices.one.marker[sapply(indices.one.marker, length) < 3]
		size.of.event <- sapply(lapply(split(x, cs), range, na.rm=TRUE), diff)
		median.size.of.event <- median(size.of.event)
		events[j, ] <- c(number.events,
				 median.number.markers.per.event,
				 median.size.of.event)
	}
	colnames(events) <- c("n_events", "median_n_markers", "median_size")
	rownames(events) <- sampleNames(object)
	return(events)
}

binary2RangedData <- function(x, pos, id, chr){
	if(all(x == 0, na.rm=TRUE)) return(RangedData())
	## Must do something about the NAs before applying cumsum
	## two situations:
	##  1.  NAs appearing in the middle of a segment
	##      - here it makes sense to 'fill in' the NAs.
	##  2.  NAs appearing between two segments (border)
	##      For border NA's we could have the following sequence
	##       0, 0, 0, NA, 1, 1, 1, 1, ...
	##       My approach is to remove the NA and assign the segment as follows (otherwise, it implies we have data in this region)
	##       [-, -, -] gap [-, - , -, - , ...]
	##       the middle will be a gap with no estimate of a mean and could potentially contain the true changepoint
	## test 1)
	##  x=c(0, 0, 0, NA, 0, 0, 1, 1, 1, 1)
	## pos=seq(along=xx)
	## chr=rep(1, length(xx))
	## truth:  2 segments
	##   start, end
	##1.  1, 6
	##2.  7, 10
	## test 2)
	##  x=c(0, 0, 0, NA, 1, 1, 1, 1)
	## pos=seq(along=x)
	## chr=rep(1, length(x))
	## truth:  2 segments
	##   start, end
	##1.  1, 3
	##2.  5, 8
	## Note that simply removing the NA's works for both cases.
	is.missing <- is.na(x)
	x <- x[!is.missing]
	pos <- pos[!is.missing]
	chr <- chr[!is.missing]
	boundaries <- c(0, diff(x) != 0 | diff(chr) != 0)
	cs <- cumsum(boundaries)
	tmp <- split(x, cs)
	isEvent <- which(sapply(tmp, unique) == 1)
	if(length(isEvent) < 1) return(RangedData())
	number.events <- length(isEvent)
	pos <- split(pos, cs)
	start.end <- t(sapply(pos[isEvent], range))
	number.markers <- sapply(tmp[isEvent], length)
	chr.split <- split(chr, cs)
	chr.split <- chr.split[isEvent]
	chrom <- sapply(chr.split, unique)
	id <- rep(id, length(chrom))

	ranges <- IRanges(start.end[, 1], start.end[, 2])
	tmp <- tryCatch(rd <- RangedData(ranges, chrom=chrom, number.markers=number.markers, id=id), error=function(e) NULL)
	if(is.null(tmp)){
		return(RangedData())
	}
	return(rd)
}

initializeEmptyDenovoSet <- function(object){
        trios.index <- getTriosIndex(object)
        denovoSet <- new("ExpressionSet", exprs=initializeBigMatrix(name="denovo",
                                          nr=nrow(object),
                                          nc=length(trios.index)))
        featureNames(denovoSet) <- featureNames(object)
        featureData(denovoSet) <- featureData(object)
	samples.index <- unlist(trios.index)
	offspring.index <- which(object$individ == 1 & 1:ncol(object) %in% samples.index)
	phenoData(denovoSet) <- phenoData(object)[offspring.index, ]
        sampleNames(denovoSet) <- denovoSet$Sample.Name
	return(denovoSet)
}

##callDenovo <- function(reduced_ranges){
##	stopifnot(inherits(reduced_ranges, "RangedData"))
##	data(pedigree)
##	## exclude ranges for which no event has occurred.
##	pedId <- who(reduced_ranges$id)
##	offspring_ranges <- reduced_ranges[pedId == "offspring", ]
##	denovo <- rep(NA, nrow(offspring_ranges))
##	uid <- substr(unique(as.character(offspring_ranges$id)), 1, 8)
##	stopifnot(all(substr(uid, 1,5) %in% trio.family))
##	rr <- IRanges(start(reduced_ranges), end(reduced_ranges))
##	rr.ids <- substr(as.character(reduced_ranges$id), 1, 8)
##	orr.ids <- substr(as.character(offspring_ranges$id), 1, 8)
##	ped.ids <- substr(pedigree$Sample.Name, 1, 8)
##	rr.chr <- reduced_ranges$chrom
##	for(i in seq_along(uid)){
##		##if(i %% 100 == 0) cat(i, " ")
##		cat(i, " ")
##		## construct the query seg.  Ignore unaffected offspring for now.
##		oid <- uid[i]
##		index1 <- which(rr.ids == oid)
##		query <- rr[index1, ]
##		parents <- as.character(pedigree[match(oid, ped.ids), 4:5])
##		index2 <- which(rr.ids %in% parents)
##		subject <- rr[index2, ]
##
##		query.chrom <- rr.chr[index1]
##		subject.chrom <- rr.chr[index2]
##
##		uchrom <- unique(query.chrom)
##
##		olaps <- vector("list", length(uchrom))
##		for(j in seq_along(uchrom)){
##			CHR <- uchrom[j]
##			qq <- query[query.chrom == CHR, ]
##			ss <- subject[subject.chrom == CHR, ]
##			if(length(ss) < 1) {
##				olaps[[j]] <- rep(0, length(qq))
##				next()
##			}
##			olaps[[j]] <- countOverlaps(qq, ss)
##		}
##		olaps <- unlist(olaps)
##		index <- which(orr.ids == oid)
##		stopifnot(length(index) == length(olaps))
##		## if > 0, one of the parents had a deletion in this region and it is not a denovo event
##		denovo[index] <- ifelse(olaps > 0, 0, 1)
##	}
##	return(denovo)
##}
##

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
