#' Export Integrated Network To Cytoscape
#' 
#' Attention: Cytoscape Software must be open.
#' Due to poor communication between R and Cytoscape, the exports are usually quite slow. To fasten things up, after running the function, wait until the "create view" button appears in Cytoscape.
#' 	Press it and wait for the network to render. Then INTERRUPT the function in R (by pressing the little stop sign above the console).
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

