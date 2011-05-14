assayDataDim <- function(object) {
	nms <- if (Biobase:::assayDataStorageMode(object) == "list") names(object) else ls(object)
	if ( length( nms ) == 0 ) return( NA )
	d <- dim( object[[ nms[[1]] ]])
	names(d) <- c( "Features", "Trios", "Family members (Father, Mother, Offspring)")
	d
}
