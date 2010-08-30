THISPKG <- "Beaty"
.onAttach <- function(libname, pkgname) {
  message("Welcome to Beaty version ", packageDescription(THISPKG, field="Version"))
}
