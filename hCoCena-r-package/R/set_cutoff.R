#' Set Cutoff
#' 
#' Set the correlation cutoffs for the different layers.
#' @param cutoff_vector A vector of cutoff values. The vector's length must be equal to the number of datasets. The order in which the cutoffs are set must correspond to the order in which the layers were declared in define_layers().
#' @export

set_cutoff <- function(cutoff_vector = base::c()){
	hcobject[["cutoff_vec"]] <<- cutoff_vector
}