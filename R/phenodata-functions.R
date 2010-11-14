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
	pData(crlmmSetList[[1]])[, "familyId"] <- as.character(pData(crlmmSetList[[1]])[, "familyId"])
	pData(crlmmSetList[[1]])[, "familyMember"] <- as.character(pData(crlmmSetList[[1]])[, "familyMember"])
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

getFamilyId <- function(object) substr(object$Sample.Name, 1, 5)
getIndividualId <- function(object) substr(object$Sample.Name, 7, 8)
getIndividualId2 <- function(object) substr(object$Sample.Name, 1, 8)

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
