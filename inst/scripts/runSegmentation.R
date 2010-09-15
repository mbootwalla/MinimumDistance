library(Beaty)
library(mybase)
library(DNAcopy)
library(crlmm)
outdir <- "/amber1/scratch/rscharpf/beaty_beadStudio"
if(!file.exists(outdir)) dir.create(outdir)
ldPath(outdir)
lrSet <- checkExists("lrSet", .path=outdir, .FUN=getLrSet)
lrSet$pedId <- who(sampleNames(lrSet))
if(!exists("batch")) stop("batch variable should be defined in the submit_logR_cbs.R script")
if(!exists("NN")) stop("batch size variable (NN) should be specified in submitter script")

sample.index <- splitIndicesByLength(1:7599, NN)[[batch]]
parent.index <- sample.index[lrSet$pedId[sample.index]=="father" | lrSet$pedId[sample.index]=="mother"]
offspring.index <- sample.index[lrSet$pedId[sample.index]=="offspring"]

## parents
marker.index <- which(chromosome(lrSet) == CHR & !duplicated(position(lrSet)))
open(logR(lrSet))
CNA.object <- CNA(genomdat=as.matrix(logR(lrSet)[marker.index, parent.index]),
		  chrom=chromosome(lrSet)[marker.index],
		  maploc=position(lrSet)[marker.index],
		  data.type="logratio",
		  sampleid=sampleNames(lrSet)[parent.index])
smu.object <- smooth.CNA(CNA.object)
tmp <- segment(smu.object, verbose=0, alpha=0.1)
cbs.segs1 <- cbind(tmp@output, tmp@segRows)
##cbs.segs1 <- print(tmp, showSegRows=TRUE, fullOutput=TRUE)

## offspring
CNA.object <- CNA(genomdat=as.matrix(logR(lrSet)[marker.index, offspring.index]),
		  chrom=chromosome(lrSet)[marker.index],
		  maploc=position(lrSet)[marker.index],
		  data.type="logratio",
		  sampleid=sampleNames(lrSet)[offspring.index])
smu.object <- smooth.CNA(CNA.object)
tmp <- segment(smu.object, verbose=0, alpha=0.01)
cbs.segs2 <- cbind(tmp@output, tmp@segRows)
##cbs.segs2 <- print(tmp, showSegRows=TRUE, fullOutput=TRUE)

cbs.segs <- rbind(cbs.segs1, cbs.segs2)
save(cbs.segs, file=file.path(outdir, paste("cbs.segs_chr", CHR, "_batch", batch, ".rda", sep="")))
q("no")
