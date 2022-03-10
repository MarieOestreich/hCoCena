#' Network integration
#' 
#' The previously constructed layer-specific networks are being integrated. 
#' @param mode A string, either "u", if network integration is to be done by union (default), or "i", if integration is to be done by intersection.
#'  For details please refer to the information pages provided in the repository's Wiki.
#' @param multi_edges One of “min”, “mean” or “max” resulting in the simplification of the multigraph by using the minimum, the mean or the maximum edge weight among the multiple edges, respectively.
#'  Multiple edges occur when the edge was present in more than one datset.
#' @param GFC_when_missing The value to substitute missing data in the case where some genes were not measured in all but only some of the datasets.
#' @param with Either an integer giving the number of the dataset to be used as reference (e.g., 1) or the name given to the layer. 
#'  Can be ignored when integration is done by union.
#' @export

build_integrated_network <- function(mode = "u", 
                                      with = NULL,
                                      multi_edges = "min",        
                                      GFC_when_missing = -hcobject[["global_settings"]][["range_GFC"]]){

  if(mode == "u"){
    message("Intergrating network based on union.")
    get_union()
  }else if(mode == "i"){
    message("Intergrating network based on intersection.")
    if(is.null(with)){
      stop("The 'with' parameter must be specified.")
    }
    get_intersection(with = with)
  }else{
    stop("No valid choice of 'mode'. Must be either 'u' for integration by union or 'i' for integration by intersection.")
  }
  merged_net <- igraph::graph_from_data_frame(hcobject[["integrated_output"]][["combined_edgelist"]], directed=FALSE)
  merged_net <- igraph::simplify(merged_net, edge.attr.comb=list(weight=multi_edges, "ignore"))
  new_edgelist <- base::cbind( igraph::get.edgelist(merged_net) , base::round(igraph::E(merged_net)$weight, 7)) %>% base::as.data.frame()
  base::colnames(new_edgelist) <- base::colnames(hcobject[["integrated_output"]][["combined_edgelist"]])
  hcobject[["integrated_output"]][["combined_edgelist"]] <<- new_edgelist
  hcobject[["integrated_output"]][["merged_net"]] <<- merged_net
  
  hcobject[["integrated_output"]][["GFC_all_layers"]] <<- merge_GFCs(GFC_when_missing = GFC_when_missing)
  

}



