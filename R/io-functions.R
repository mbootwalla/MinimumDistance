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
