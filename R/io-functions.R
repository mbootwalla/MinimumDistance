readPennCnv <- function(penndir="/thumper/ctsa/beaty/holger/penncnv/jointDat",
			offspring.only=TRUE){
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
		gr <- GRanges(seqnames=chr,
			      ranges=IRanges(penn.joint$StartPosition,
			      penn.joint$EndPosition),
			      nmarkers=penn.joint$NumberSNPs,
			      id=penn.joint$ID,
			      triostate=penn.joint$TrioState,
			      pedId=penn.joint$FamilyMember)
		rdList[[i]] <- gr
	}
	penn.joint <- do.call("c", rdList)
	penn.joint <- as(penn.joint, "RangedData")
	if(offspring.only){
		message("Returning only the ranges for the offspring")
		penn.joint <- penn.joint[penn.joint$pedId=="offspring", ]
		chrom <- as.character(space(penn.joint))
		chrom <- substr(chrom, 4, nchar(chrom))
		penn.joint$chrom <- chromosome2integer(chrom)
	}
##	pennJoint <- do.call("c", rdList)
##	return(pennJoint)
	return(penn.joint)
}
