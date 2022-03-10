
#' Plot of Cutoff Statistics
#' 
#' Plots the R-squared value, number of edges, number of genes and number of networks for different cut-offs.
#' @param interactive Boolean. If TRUE (default) plot is created using plotly including and interactive cutoff slider. If FALSE, plot is created as a static ggplot.
#' @param hline A list with four slots ("R.squared", "no_edges", "no_nodes", "no_networks") each of which can be set either to NULL (default) or a number to introduce a horizontal line for orientation at that value in the respective plot. Only used in the non-interactive plot.
#' @export

plot_cutoffs <- function(interactive = T, 
                         hline = list("R.squared" = NULL, 
                                      "no_edges" = NULL, 
                                      "no_nodes" = NULL, 
                                      "no_networks" = NULL)){
	
              
	for(x in 1:base::length(hcobject[["layers"]])){
		cutoff_df <- hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["cutoff_calc_out"]][["cutoff_stats_concise"]]
		if(interactive == T){
		  p <- plot_cutoffs_internal_interactive(cutoff_stats = cutoff_df,
		                                         x = x)
		}else{
		  p <- plot_cutoffs_internal_static(cutoff_stats = cutoff_df, 
		                                    hline = hline, 
		                                    x = x)
		}

	print(p)
	hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["cutoff_plot"]] <<- p

	}
}	