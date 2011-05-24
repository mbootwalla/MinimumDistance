setMethod("family", "RangedData", function(object) substr(object$id, 1, 5))


