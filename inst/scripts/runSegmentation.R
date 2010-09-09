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
marker.index <- which(chromosome(lrSet) == 1 & !duplicated(position(lrSet)))
CNA.object <- CNA(genomdat=as.matrix(logR(lrSet)[marker.index, sample.index]),
		  chrom=chromosome(lrSet)[marker.index],
		  maploc=position(lrSet)[marker.index],
		  data.type="logratio",
		  sampleid=sampleNames(lrSet)[sample.index])
smu.object <- smooth.CNA(CNA.object)
tmp <- segment(smu.object)
cbs.segs <- print(tmp, showSegRows=TRUE)
save(cbs.segs, file=file.path(outdir, paste("cbs.segs_", batch, ".rda", sep="")))
##rD <- cbs(lrSet, sample.index=sample.index)
##save(rD, file=file.path(outdir, paste("rD_", batch, ".rda", sep="")))
q("no")
