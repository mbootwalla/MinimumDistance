##not that I know. However, if you are just interested in the logR and the
##BAF values they are stored in the files in
##/thumper/ctsa/beaty/holger/txtfiles (one file per sample), and the file
##names are the sample names (so <familyID>_<personid>@<someIlluminaNumber>)
##as stored in the third column of the Beaty_... files. The genotypes are
##stored in /thumper/ctsa/beaty/holger/genotypes.
library(Beaty)
library(ff)
library(crlmm)
path <- "/thumper/ctsa/beaty/holger/txtfiles"
fnames <- list.files(path, full.names=TRUE)
sns <- as.character(sapply(basename(fnames), function(x) strsplit(x, ".txt")[[1]][[1]]))
data(samplesheet, package="Beaty")
tmp <- match(sns, samplesheet$Sample.Name)
outdir <- "/amber1/scratch/rscharpf/beaty_beadStudio"
if(!file.exists(outdir)) dir.create(outdir)
ldPath(outdir)
dat <- read.delim(fnames[1])
ocSamples(100)

if(FALSE){
	baf <- initializeBigMatrix("baf", nrow(dat), length(sns), vmode="double")
	save(baf, file=file.path(outdir, "baf.rda"))
	logR <- initializeBigMatrix("logR", nrow(dat), length(sns), vmode="double")
	save(logR, file=file.path(outdir, "logR.rda"))

	tmp <- new("MultiSet", logR=logR, annotation="human610quadv1b")
	featureNames(tmp) <- as.character(dat[, 1])
	featureData(tmp) <- addFeatureAnnotation(tmp)
	ix <- order(chromosome(tmp), position(tmp))
	save(ix, file=file.path(outdir, "ix.rda"))

	nr <- nrow(tmp); nc <- ncol(tmp)
	bsSet <- new("MultiSet",
		     logR=logR,
		     BAF=baf,
		     segmentMeans=initializeBigMatrix(name="segMeans", nr, nc, vmode="double"),
		     call=initializeBigMatrix(name="beadStudioCalls", nr, nc, vmode="integer"),
		     annotation=annotation(tmp))
	save(bsSet, file=file.path(outdir, "bsSet.rda"))
	featureNames(bsSet) <- featureNames(tmp)[ix]
	featureData(bsSet) <- addFeatureAnnotation(bsSet)
	sampleNames(bsSet) <- sns
	data(samplesheet, package="Beaty")
	samplesheet <- samplesheet[match(sampleNames(bsSet), samplesheet$Sample.Name), ]
	stopifnot(identical(samplesheet$Sample.Name, sampleNames(bsSet)))	
	pData(bsSet) <- samplesheet
	save(bsSet, file=file.path(outdir, "bsSet.rda"))
} else{
	load(file.path(outdir, "ix.rda"))
	load(file.path(outdir, "bsSet.rda"))
	logR <- assayData(bsSet)[["logR"]]
	baf <- assayData(bsSet)[["BAF"]]
	open(logR)
	open(baf)
}

batches <- splitIndicesByLength(seq(along=sns), ocSamples())
for(j in batches){
	cat(".")
	tmpBaf <- matrix(NA, nrow(dat), length(j))
	tmpLogr <- matrix(NA, nrow(dat), length(j))
	for(k in seq(along=j)){
		##cat(".")
		if(k %% 100 == 0) cat(".")
		if(!file.exists(fnames[j[k]])) next()
		dat <- read.delim(fnames[j[k]])
		##order baf and logr values by chromosome, physical position
		tmpBaf[, k] <- dat[ix, 3]
		tmpLogr[, k] <- dat[ix, 2]
	}
	baf[, j] <- tmpBaf
	logR[, j] <- tmpLogr
	rm(tmpBaf)
	rm(tmpLogr)
	gc()
}
assayDataElementReplace(bsSet, "logR", logR)
assayDataElementReplace(bsSet, "BAF", baf)
save(bsSet, file=file.path(outdir, "bsSet.rda"))
close(logR)
close(baf)


