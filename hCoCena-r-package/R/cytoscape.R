#' Export Integrated Network to local R session
#' Due to a lacking possibility of communication with Cytoscape from within Docker container, you need to export all necessary information for import into a local R session.
#' @param file Path to the folder where the network information should be saved. Default path is set to the save folder defined in the global settings.
#' @export

export_to_local_folder <- function(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]])){
  
  network_df <- igraph::as_data_frame(network_filt())
  network_df$weight <- NULL
  
  readr::write_delim(network_df, file = paste0(file, "/network.txt"), 
                     delim = "\t", 
                     col_names = T)
}



#' Export Integrated Network To Cytoscape
#' 
#' Attention: Cytoscape Software must be open.
#' Due to difficulties in the communication between R/RCy3 and Cytoscape, you need to manually stop the function in R as soon as the table of nodes and edges appears in Cytoscape.(by pressing the little stop sign above the console).
#' @param name A string. The name given to the graph in Cytoscape. Default is "my igraph".
#' @export

export_to_cytoscape <- function(name = "my igraph", docker_container = FALSE){
    
    RCy3::createNetworkFromIgraph(network_filt(), name)

}


#' Import Layout From Cytoscape
#' 
#' Imports the layout of a network currently open in Cytoscape.
#' @export

import_layout_from_cytoscape <- function(){
  
    l <- RCy3::getNodePosition() %>% as.matrix()
    rn <- rownames(l)
    l <- apply(l, 2, as.numeric)
    rownames(l) <- rn
    
    hcobject[["integrated_output"]][["cluster_calc"]][["layout"]] <<- l

}



#' Import network layout from file to your Docker container
#' Imports the layout generated within a local R session and Cytoscape back into the Docker container
#' @param file Exact path and file name containing the Cytoscape layout. Default path is set to the save folder defined in the global settings.
#' @export

import_layout_from_local_folder <- function(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/network_layout.csv")){
  
  l <- utils::read.csv(file = file, 
                       row.names = 1)
  
  rn <- rownames(l)
  l <- apply(l, 2, as.numeric)
  rownames(l) <- rn
  
  hcobject[["integrated_output"]][["cluster_calc"]][["layout"]] <<- l
  
}


