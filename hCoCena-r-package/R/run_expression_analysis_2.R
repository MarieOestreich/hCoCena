#' Run Second Part Of The Gene Expression Analysis
#' 
#' This function plots a heatmap for the network genes in each data layer and computes the Group-Fold-Changes for all genes per layer.
#' @param grouping_v A string giving a column name present in all annotation files, if this variable shall be used for grouping the samles isntead of the variable of interest. Default is NULL.
#' @param plot_HM A Boolean. Whether or not to plot the heatmap (for networks with many genes this may be very demanding for your computer if you are running the analysis locally). Default is TRUE.
#' @param method The method used for clustering the heatmap in the pheatmap function. Default is "complete".
#' @param additional_anno A list, with one slot per data set. A slot contains a vector of column names from that data set’s annotation file that you wish to annotate with. 
#' 	If for some of the data sets you don’t wish any further annotation, you can set the corresponding list slot to NULL. Default is NULL.
#' @export

run_expression_analysis_2 <- function(grouping_v = NULL, plot_HM = T, method = "complete", additional_anno = NULL){
  
  for(x in 1:base::length(hcobject[["layers"]])){
    
    run_expression_analysis_2_body(x = x, 
    							   grouping_v = grouping_v, 
							       plot_HM = plot_HM, 
							       method = method, 
							       additional_anno = additional_anno[[x]])
  }
}


