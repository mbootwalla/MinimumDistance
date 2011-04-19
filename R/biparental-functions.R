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
