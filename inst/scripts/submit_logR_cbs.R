library(oligoClasses)
NN <- 600
sampleBatch <- splitIndicesByLength(1:7599, NN)
for(KK in 1:length(sampleBatch)){
	sink("temp")
	cat("runSegmentation <- TRUE \n")
	cat("sourceSubmitter <- FALSE \n")
	cat("batch <- ", KK, "\n")
	cat("NN <- ", NN, "\n")
	sink()
	fn <- paste("tmp_", KK, ".R", sep="")
	if(file.exists(fn)) unlink(fn)
	system(paste("cat temp runSegmentation.R >", fn))
	system(paste("clusterMcmc 5G", fn))
	Sys.sleep(45)
}
