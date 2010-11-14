beadstudiodir <- function(string){
	if (missing(string)){
		return(getOption("beadstudiodir"))
	}else{
		options(beadstudiodir=string)
		invisible(TRUE)
	}
}

crlmmdir <- function(string){
	if (missing(string)){
		return(getOption("crlmmdir"))
	}else{
		options(crlmmdir=string)
		invisible(TRUE)
	}
}

beatyOptions <- function(crlmm.path="/amber1/archive/ctsa/beaty_crlmm/devel",
			 bs.path="/amber1/archive/ctsa/beaty_crlmm/beaty_beadStudio"){
	beadstudiodir(bs.path)
	crlmmdir(crlmm.path)
	TRUE
}
