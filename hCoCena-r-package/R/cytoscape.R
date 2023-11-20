#' Export Integrated Network To Cytoscape
#' 
#' Attention: Cytoscape Software must be open.
#' Due to difficulties in the communication between R/RCy3 and Cytoscape, you need to manually stop the function in R as soon as the table of nodes and edges appears in Cytoscape.(by pressing the little stop sign above the console).
#' @param name A string. The name given to the graph in Cytoscape. Default is "my igraph".
#' @export

export_to_cytoscape <- function(name = "my igraph"){

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

