THISPKG <- "Beaty"
.beatyEnv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname) {
	require("methods")
}

.onAttach <- function(libname, pkgname) {
	version <- packageDescription("Beaty", field="Version")
	message(getBar())
	message("Welcome to Beaty version ", version)

	options(prompt="R> ", continue=" ", width=70)##device=pdf
	##library(Study)
##	library(oligoClasses)
##	library(BeatyExperimentData)
#3	library(SNPchip)
#3	library(DNAcopy)
#3	library(IRanges)
#3	library(mybase)
#3	library(ff)
#3	library(crlmm)
	##library(GenomicRanges)
	##library(Study)
##	library(lattice)
##	library(grid)
##	library(locuszoom)
##	library(RColorBrewer)
	ocSamples(100)
	beatyOptions()
}

.onUnload <- function(libpath){
	library.dynam.unload(THISPKG, libpath)
}


