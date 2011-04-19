setReplaceMethod("copyNumber", signature(object="CopyNumberSet",
					 value="ff_matrix"), function(object, value){
						 assayDataElementReplace(object, "copyNumber", value)
					 })

