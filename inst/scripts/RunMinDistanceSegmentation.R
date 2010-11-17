library(Beaty)
library(BeatyExperimentData)
outdir <- beadstudiodir()
ldPath(outdir)
data(minDistanceSet)
data(bsSet)
featureData(minDistanceSet) <- featureData(bsSet)
rm(bsSet); gc()
j <- splitIndicesByLength(1:ncol(minDistanceSet), BATCHSIZE)[[BATCH]]
marker.index <- which(chromosome(minDistanceSet) == CHR)
marker.index <- marker.index[!duplicated(position(minDistanceSet)[marker.index])]
stopifnot(all(diff(order(chromosome(minDistanceSet)[marker.index], position(minDistanceSet)[marker.index])) >= 0))
invisible(open(copyNumber(minDistanceSet)))
CNA.object <- CNA(genomdat=as.matrix(copyNumber(minDistanceSet)[marker.index, j]),
		  chrom=chromosome(minDistanceSet)[marker.index],
		  maploc=position(minDistanceSet)[marker.index],
		  data.type="logratio",
		  sampleid=sampleNames(minDistanceSet)[j])
close(copyNumber(minDistanceSet))
smu.object <- smooth.CNA(CNA.object)
tmp <- segment(smu.object, verbose=0)
md.segs <- cbind(tmp$output, tmp$segRows)
fname <- paste("md.segs_chr", CHR, "_batch", BATCH, ".rda", sep="")
save(md.segs, file=file.path(ldPath(), fname))
q("no")
