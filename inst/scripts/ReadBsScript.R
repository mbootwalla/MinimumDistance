stopifnot(exists("J"))

library(oligoClasses)
ocSamples(500)
library(mybase)
outdir <- "/amber1/scratch/rscharpf/beaty_beadStudio"
load(file.path(outdir, "bsSet.rda"))
sns <- sampleNames(bsSet)

path <- "/thumper/ctsa/beaty/holger/txtfiles"
fnames <- list.files(path, full.names=TRUE)

fnames.short <- as.character(sapply(basename(fnames), function(x) strsplit(x, ".txt")[[1]][[1]]))
table(fnames.short %in% sns)
tmp <- fnames[fnames.short %in% sns]
fnames.short <- as.character(sapply(basename(tmp), function(x) strsplit(x, ".txt")[[1]][[1]]))
index <- match(sns, fnames.short)
stopifnot(all.equal(fnames.short[index], sns))
fnames <- tmp[index]

dat <- read.delim(fnames[1], colClasses=c("character", "numeric", "numeric"))
match.index <- match(featureNames(bsSet), dat$Name)
j <- splitIndicesByLength(seq(along=fnames), ocSamples())[[J]]

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
bsSet$MAD[j] <- mads

close(baf(bsSet))
close(logR(bsSet))
q("no")
