#!/usr/bin/env Rscript

## [[file:~/projects/scripts.org::*BeadStudio%20CBS][block-6]]

library(mybase)
 library(Beaty)
 library(ff)
 library(crlmm)
 library(IRanges)
 library(SNPchip)
 library(DNAcopy)
 data(samplesheet, package="Beaty")
 outdir <- "/amber1/scratch/rscharpf/beaty_beadStudio"
 if(!file.exists(outdir)) dir.create(outdir)
 ldPath(outdir)
 if(FALSE){
        load(file.path(outdir, "bsSet.rda"))
        lrSet <- as(bsSet, "LogRatioSet")
        save(lrSet, file=file.path(outdir, "lrSet.rda"))
  } else load(file.path(outdir, "lrSet.rda"))
  if(!exists(runSegmentation))  runSegmentation <- FALSE
  if(!exists(sourceSubmitter))  sourceSubmitter <- FALSE
  if(sourceSubmitter){
          ##the only way the submitter is
          warning("sourcing script to segment the logR values")
          source("submit_logR_cbs.R")
  }  
  if(runSegmentation){
        if(!exists("batch")) stop("batch variable should be defined in the submit_logR_cbs.R script")
        if(!exists("NN")) stop("batch size variable (NN) should be specified in submitter script")      
        sample.index <- splitIndicesByLength(1:7599, NN)[[batch]]                                 
        rD <- cbs(lrSet, sample.index=sample.index)
        save(rD, file=file.path(outdir, paste("rD_", batch, ".rda", sep="")))
        q("no")
  }
  ## Run cbs on logR values
  ##library(mybase, lib="~/Rlibs/devel2")
  library(Beaty, lib="~/Rlibs/devel2")
  library(ff)
  library(oligoClasses)
  library(crlmm)
  library(DNAcopy)
  library(IRanges)
  data(samplesheet, package="Beaty")
  outdir <- "/amber1/scratch/rscharpf/beaty_beadStudio"
  ldPath(outdir)
  load(file.path(outdir, "bsSet.rda"))
  source("~/projects/mybase/R/functions.R")
  trace(cbs, browser)
  cbsResults <- cbs(bsSet)
  
  logR <- assayData(bsSet)[["logR"]]
  mads <- rep(NA, 1:ncol(bsSet))
  for(j in seq_along(mads)) mads[j] <- mad(logR[, j], na.rm=TRUE)
  pData(bsSet)$MAD <- mads
  pData(bsSet) <- getFamilyInfo(pData(bsSet))
  save(bsSet, file=file.path(outdir, "bsSet.rda"))
  
  chr <- chromosome(bsSet)
  pos <- position(bsSet)
  mapped <- !is.na(pos)
  Y <- chromosome(bsSet) == 24
  chr <- chr[mapped & !Y]
  pos <- pos[mapped & !Y]
  SNP <- isSnp(bsSet)[mapped & !Y]
  segmentMeans <- assayData(bsSet)[["segmentMeans"]]
  open(segmentMeans)
  ##cT <- cT[mapped & !Y]
  res <- list()
  ## NN, KK defined in submit_logR_cbs.R 
  ##sampleBatch <- splitIndicesByLength(1:ncol(bsSet), NN)
  ##batch <- sampleBatch[[KK]]
  J <- which(segmentMeans[1, ] == 0)
  JJ <- splitIndicesByLength(J, NN)
  batch <- JJ[[KK]]
  for(j in batch){
          cat(j, " ")
          lR <- logR[mapped & !Y, j]
          ##smoothed logR vals
          obj <- CNA(lR, chrom=chr, maploc=pos)
          obj <- smooth.CNA(obj)
          system.time(test <- segment(obj))
          tmp <- rep(test$output$seg.mean, test$output$num.mark)
          index <- which(mapped & !Y)[!is.na(lR)]
          segmentMeans[index, j] <- tmp
          res[[j]] <- test$output
  }
  save(res, file=file.path(outdir, "res.rda"))
  q("no")
## block-6 ends here
