#' Check Transcription Factor
#' 
#' This function leverages the information collected with TF_enrich_network()
#' 	to allow the user to query specific transcription factors of interest and see how their top targets are spread across modules.  
#' 	The goal is to uncover potential co-regulations between clusters.
#' @param TF A string giving the name of the transcription factor to be queried.
#' @export

check_tf <- function(TF){


	# get edgelist from integrated network:
	edgelist <- hcobject[["integrated_output"]][["combined_edgelist"]][, base::c(1,2)]
	edgelist[] <- base::lapply(edgelist, as.character)
	edgelist$merged <- base::paste0(edgelist$V1, edgelist$V2)
	edgelist$merged2 <- base::paste0(edgelist$V2, edgelist$V1)

	gtc <- GeneToCluster()
	base::colnames(gtc) <- base::c("gene", "cluster")

	# the targets of the transcription factor in question:
	targets <- hcobject[["integrated_output"]][["enrichall"]][[TF]][["targets"]]

	# create edgelist from TF to it's targets, removing self edges:
	edges <- base::data.frame(from = base::rep(TF, base::length(targets)), to = targets) %>% 
		base::unique()
	edges[] <- base::lapply(edges, as.character)
	edges <- edges[!edges$to == TF,]
	merged <- base::paste0(edges$from, edges$to)
	# if edge exists in network, colour = black, otherwise color = grey:
	edges$color <- base::lapply(merged, function(x){
	if(x %in% edgelist$merged | x %in% edgelist$merged2){
	  "black"
	}else{
	  "grey"
	}
	})%>% base::unlist()

	# get nodes (TF and targets that are connected in the network) including their cluster colour:
	nodes <- base::data.frame(name = base::unique(base::c(edges$from, edges$to)))
	nodes <- base::merge(nodes, gtc, by.x = "name", by.y = "gene") 
	base::colnames(nodes) <- base::c("name", "color")
	nodes <- nodes[base::order(nodes$color),]

	# create star plot with TF at the center:
	g <- igraph::graph_from_data_frame(d = edges, vertices = nodes$name)
	igraph::V(g)$color <- base::as.character(nodes$color)
	l <- igraph::layout.star(g, center = igraph::V(g)[TF])
	#plot to PDF:
	Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/TF_", TF, "_starplot.pdf"), width = 10,
	                            height = 10)
	igraph::plot.igraph(g, layout = l, edge.arrow.size = 0.5, vertex.label.color = "black", edge.color = edges$color,
	                  vertex.label.cex = 0.7, vertex.label.font = 2, edge.width = 2, 
	                  vertex.frame.color = base::as.character(nodes$color))
	grDevices::dev.off()

	# plot to markdown:
	igraph::plot.igraph(g, layout = l, edge.arrow.size = 0.5, vertex.label.color = "black", edge.color = edges$color,
	                  vertex.label.cex = 0.7, vertex.label.font = 2, edge.width = 2, 
	                  vertex.frame.color = base::as.character(nodes$color))
  
}