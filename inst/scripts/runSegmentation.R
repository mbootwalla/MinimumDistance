library(mybase)
library(DNAcopy)
library(crlmm)
outdir <- "/amber1/scratch/rscharpf/beaty_beadStudio"
if(!file.exists(outdir)) dir.create(outdir)
ldPath(outdir)
lrSet <- checkExists(name="lrSet", path=outdir, FUN=getLrSet)
lrSet$pedId <- who(sampleNames(lrSet))
if(!exists("batch")) stop("batch variable should be defined in the submit_logR_cbs.R script")
if(!exists("NN")) stop("batch size variable (NN) should be specified in submitter script")      
sample.index <- splitIndicesByLength(1:7599, NN)[[batch]]                                 
rD <- cbs(lrSet, sample.index=sample.index)
save(rD, file=file.path(outdir, paste("rD_", batch, ".rda", sep="")))
q("no")
