#' Export Integrated Network to local R session
#' Due to a lacking possibility of communication with Cytoscape from within Docker container, you need to export all necessary information for import into a local R session.
#' @export

export_to_localR <- function(){
  
  gtc <- GeneToCluster()
  network <- hcobject[["integrated_output"]][["merged_net"]]
  
  base::save(gtc, network, file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/export_to_cytoscape.RData"))
}



#' Export Integrated Network To Cytoscape
#' 
#' Attention: Cytoscape Software must be open.
#' Due to difficulties in the communication between R/RCy3 and Cytoscape, you need to manually stop the function in R as soon as the table of nodes and edges appears in Cytoscape.(by pressing the little stop sign above the console).
#' @param name A string. The name given to the graph in Cytoscape. Default is "my igraph".
#' @param docker_container A boolean. If you exported the imported object from within a docker container, set to TRUE. FALSE by default. 
#' @export

export_to_cytoscape <- function(name = "my igraph", docker_container = FALSE){
  
  if(docker_container){
    
    del_v <- igraph::V(network)$name[!igraph::V(network)$name %in% gtc$gene]
    network_filt <- igraph::delete.vertices(network, del_v)

    RCy3::createNetworkFromIgraph(network_filt(), name)
    
    
  } else {
    
    RCy3::createNetworkFromIgraph(network_filt(), name)

  }
}


#' Import Layout From Cytoscape
#' 
#' Imports the layout of a network currently open in Cytoscape.
#' @param docker_container A boolean. If you exported the imported object from within a docker container, set to TRUE. FALSE by default. 
#' @export

import_layout_from_cytoscape <- function(docker_container = FALSE, path){
  
  if(docker_container){
    
    
  } else {
    
    l <- RCy3::getNodePosition() %>% as.matrix()
    rn <- rownames(l)
    l <- apply(l, 2, as.numeric)
    rownames(l) <- rn
    
    base::saveRDS(l, file = paste0(path, "/cytoscape_layout.rds"))
    
  }

	l <- RCy3::getNodePosition() %>% as.matrix()
	rn <- rownames(l)
	l <- apply(l, 2, as.numeric)
	rownames(l) <- rn

	hcobject[["integrated_output"]][["cluster_calc"]][["layout"]] <<- l

}



#' Import network layout from file to your Docker container
#' Imports the layout generated within a local R session and Cytoscape back into the Docker container

import_from_localR <- function(){
  
  l <- base::loadRDS(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/cytoscape_layout.rds"))
  hcobject[["integrated_output"]][["cluster_calc"]][["layout"]] <<- l
  
}


