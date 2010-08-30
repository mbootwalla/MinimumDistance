##setClass("TrioSet", contains="SnpCallSetPlus", representation(isBiparental="logical"))
##setGeneric("isBiparental", function(object) standardGeneric("isBiparental"))
####setMethod("isBiparental", "TrioSet", function(object) object@isBiparental)
##setMethod("initialize", "TrioSet", function(.Object,
##					    isBiparental, ...){
##	.Object <- callNextMethod(.Object, ...)
##	.Object@isBiparental <- isBiparental
##})


