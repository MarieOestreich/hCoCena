#' Defines the Data Layers
#' 
#' The function collects descriptive names of the datasets as well as the names of the count and annotation files for each dataset.
#' @param data_sets A named list. Each list entry represents a dataset and must have a name and a value. The name should describe the dataset and will be used for plot labels etc., while the entry should be a vector or exactly two strings.
#' 	The first string should be the name of the dataset's count file (incl. file ending), the second the name of the dataset's annotation file (incl. file ending).
#' 	Alternatively to providing the count and annotation data as files, they can also be provided as existing dataframe objects. 
#' 	Give the object names as strings instead of the file names. The count object must be a data frame with sample names as columns and gene names as rows. 
#' 	There must be no additional columns other than those representing the counts per sample. 
#' 	The annotation object must also be a data frame, where row names are sample names that match the column names of the count object, and column names are information categories.
#' @export
#' @examples
#' define_layers(data_sets = list(rhinovirusSet = c("rhinovirusSetCount.txt", "rhinovirusSetAnno.txt"),
#' 								  influenzaSet = c("influenzaSetCount.txt", "influenzaSetAnno.txt")))
#' or from data frames:
#' define_layers(data_sets = list(rhinovirusSet = c("rhinovirusSetCountDf", "rhinovirusSetAnnoDf"),
#' 								  influenzaSet = c("influenzaSetCountDf", "influenzaSetAnnoDf")))

define_layers <- function(data_sets = list()){

	hcobject[["layers"]] <<- list()
	for(setnum in 1:base::length(data_sets)){
		hcobject[["layers"]][[base::paste0("set", setnum)]] <<- data_sets[[setnum]]
	}

	hcobject[["layers_names"]] <<- base::names(data_sets)
	
}