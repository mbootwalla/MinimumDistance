setMethod("chromosome", "GRanges", function(object) {
	as.integer(sapply(as.character(seqnames(object)), function(x) strsplit(x, "chr")[[1]][2]))
})

setMethod("$", "GRanges", function(x, name) {
	eval(substitute(elementMetadata(x)$NAME_ARG, list(NAME_ARG=name)))
})

setMethod("colnames", "GRanges", function(x, do.NULL=TRUE, prefix="col") {
	colnames(elementMetadata(x))
})

setReplaceMethod("$", "GRanges", function(x, name, value) {
	elementMetadata(x)[, name] = value
	x
})
