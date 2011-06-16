THISPKG <- "MinimumDistance"
.mdEnv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname) {
	require("methods")
}

.onAttach <- function(libname, pkgname) {
	version <- packageDescription("MinimumDistance", field="Version")
	message(getBar())
	message("Welcome to MinimumDistance version ", version)
	options(prompt="R> ", continue=" ", width=70)##device=pdf
	ocSamples(100)
}

.onUnload <- function(libpath){
	library.dynam.unload(THISPKG, libpath)
}


