#' Define Layer Settings
#' 
#' Receives all settings that are dataset specific and not globally valid.
#' @param top_var A vector with the length equal to the number of datasets/layers. Each entry of the vector is either "all" or an integer. 
#' 	Defines the number of most variable genes to be extracted per dataset. The order in the vector corresponds to the order in which the datasets have been declared in define_layers().
#' @param min_corr A vector of floats with the length equal to the number of datasets/layers. To construct a meaningful co-expression network for each layer, correlation cut-offs must be determined for every dataset that mark the lower boundary for the correlation 
#' 	of two genes in order for their co-expression to be represented as an edge in the network. To facilitate choosing these cut-offs, a series of parameters will be calculated for a defined number of different cut-offs. 
#' 	The number of cut-offs for which these parameters are calculated is determined by 'range_cutoff_length'. 
#' 	The range from which these possible cut-offs are taken is on the lower end restricted by 'min_corr' and on the upper end by the maximum correlation calculated between any two genes in the dataset. Default is 0.7.
#' @param range_cutoff_length A vector of integers with the length equal to the number of datasets/layers. Details see "min_corr".
#' @param print_distribution_plots A vector of Booleans with the length equal to the number of datasets/layers. Whether or not to print the degree distribution plots for all tested cut-offs to pdf files. 
#' 	The number of plots per data set will therefore be equal to the 'range_cutoff_length' parameter you have set. 
#' 	Given the potential size, this should only be set to TRUE, if 'range_cutoff_length' is small or if one wishes to thoroughly analyse how the degree distribution changes in detail for differing cut-offs.
#' 	Default is FALSE.
#' @export

set_layer_settings <- function(top_var, 
                               min_corr = 0.7, 
                               range_cutoff_length, 
                               print_distribution_plots = FALSE){

	if(base::length(print_distribution_plots) == 1 & print_distribution_plots[1] == FALSE){
		print_distribution_plots <- base::rep(FALSE, base::length(hcobject[["layers"]]))
	}
	  
	for(x in 1:base::length(hcobject[["layers"]])){

		hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]] <<- top_var[x]
		hcobject[["layer_settings"]][[base::paste0("set", x)]][["min_corr"]] <<- min_corr[x]
		hcobject[["layer_settings"]][[base::paste0("set", x)]][["range_cutoff_length"]] <<- range_cutoff_length[x]
		hcobject[["layer_settings"]][[base::paste0("set", x)]][["print_distribution_plots"]] <<- print_distribution_plots[x]

	}
  
}