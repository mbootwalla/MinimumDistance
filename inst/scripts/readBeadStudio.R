#!/usr/bin/env Rscript

## [[file:~/projects/Beaty/inst/scripts/scripts.org::*Read%20BeadStudio%20output][block-2]]

##Log R values
##Remove first 10 rows
##extract 2nd column
bsdir <- "/thumper/ctsa/beaty/illumina610/Genotype_and_Intensity_Data_Files/Beaty_610Q_genotype_files_0409"
beadStudioFiles <- list.files(bsdir, full.names=TRUE)

if(FALSE){
        If(!file.exists("~/projects/Beaty/inst/extdata/key.rda")){
                key <- matrix(NA, nrow=length(beadStudioFiles), ncol=2)
                colnames(key) <- c("beadStudioName", "Sample.ID")
                outfile <- "/thumper/ctsa/beaty/scharpf/tmp.csv"
                outfile2 <- "/thumper/ctsa/beaty/scharpf/tmp2.txt"
                outfile3 <- "/thumper/ctsa/beaty/scharpf/tmp3.txt"
                for(i in seq(along=beadStudioFiles)){
                        cat(".")
                        system(paste("sed -n -e '12p' ", beadStudioFiles[i], " > /thumper/ctsa/beaty/scharpf/tmp.csv", sep=""))
                        tmp <- as.character(read.csv(outfile, header=FALSE, as.is=TRUE)[2])
                        key[i, ] <- c(basename(beadStudioFiles[i]), tmp)
                        unlink(outfile)
                }
                save(key, file="~/projects/Beaty/inst/extdata/key.rda")
        } else {
                message("Loading key.rda...")
                load("~/projects/Beaty/inst/extdata/key.rda")
        }
        q("no")
}


if(!file.exists("/thumper/ctsa/beaty/scharpf/position.rda")){
        library(oligoClasses)
        tmp <- read.csv("/thumper/ctsa/beaty/illumina610/SNP_Summary_Files/Beaty_610Q_release_SNP Table.csv")
        position <- matrix(NA, 620901, 2)
        colnames(position) <- c("chromosome", "position")
        position[, "chromosome"] <- chromosome2integer(as.character(tmp$Chr))
        position[, "position"] <- tmp$Position
        rownames(position) <- tmp$Name
        save(position, file="/thumper/ctsa/beaty/scharpf/position.rda")
} else {
        message("loading genomic position ...")
        load("/thumper/ctsa/beaty/scharpf/position.rda")
}
if(FALSE){
        beadStudioFilesList <- split(beadStudioFiles, rep(1:10, each=length(beadStudioFiles)/10, length.out=length(beadStudioFiles)))
        ##save by chromosome
        outdir <- file.path("/thumper/ctsa/beaty/scharpf/beadstudio", k)
        if(!file.exists(outdir)) dir.create(outdir)
        beadStudioFiles <- beadStudioFilesList[[k]]
        beadStudioData <- read.csv(beadStudioFiles[1], skip=10, as.is=TRUE)
        logR.all <- list()
        BAF.all <- list()
        logR.all[[1]] <- beadStudioData[, "Log.R.Ratio"]
        BAF.all[[1]] <- beadStudioData[, "B.Allele.Freq"]
        snps <- beadStudioData[, "SNP.Name"]
        identical(snps, rownames(position))

        for(j in 2:length(beadStudioFiles)){
                ##for(j in 2:10){
                cat(".")
                beadStudioData <- read.csv(beadStudioFiles[1], skip=10, as.is=TRUE)
                logR.all[[j]] <- as.integer(beadStudioData[, "Log.R.Ratio"] * 1000)
                BAF.all[[j]] <- as.integer(beadStudioData[, "B.Allele.Freq"]*1000)
                rm(beadStudioData); gc()
        }

for(j in 1:23){
        cat("Chromosome ", j, "\n")
        isChrom <- position[, "chromosome"] == j
        tmp <- lapply(logR.all, function(x) x[isChrom])
        logR <- do.call(cbind, tmp)
        colnames(logR) <- basename(beadStudioFiles)
        rownames(logR) <- snps[isChrom]
        save(logR, file=file.path(outdir, paste("logR_", j, ".rda", sep="")))
        rm(logR, tmp); gc()
        
        tmp <- lapply(BAF.all, function(x) x[isChrom])
        BAF <- do.call(cbind, tmp)
        colnames(BAF) <- basename(beadStudioFiles)
        rownames(BAF) <- snps[isChrom]
        save(BAF, file=file.path(outdir, paste("BAF_", j, ".rda", sep="")))
        rm(BAF, tmp); gc()
}
q("no")
}



##fns <- list.files("~/projects/Beaty/inst/extdata", full.names=TRUE)
##tmp <- read.csv(fns[1])
##620,901 snps  (is this snp-only?)
## block-2 ends here
