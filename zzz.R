.onAttach <- function(libname, pkgname) {
	version <- packageDescription("MinimumDistance", field="Version")
	message("Welcome to MinimumDistance version ", version)
}
