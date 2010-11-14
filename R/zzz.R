THISPKG <- "Beaty"
.beatyEnv <- new.env(parent=emptyenv())

.onAttach <- function(libname, pkgname) {
	version <- packageDescription("Beaty", field="Version")
	message(getBar())
	message("Welcome to Beaty version ", version)

	options(prompt="R> ", continue=" ", width=70)##device=pdf
	library(Study)
	library(oligoClasses)
	library(BeatyExperimentData)
	library(SNPchip)
	library(DNAcopy)
	library(IRanges)
	library(mybase)
	library(ff)
	library(crlmm)
	library(GenomicRanges)
	library(Study)
	ocSamples(100)
	beatyOptions()
}


