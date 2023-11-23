#' Plot GFC Networks
#' 
#' For visualization of the GFC for every gene under the different observed groups, the network can be additionally replotted once for every group, 
#' 	with nodes being coloured according to their GFC value. 
#' 	This provides a more detailed resolution of the information acquired from the module heatmap.
#' @export


plot_GFC_network <- function(){

  GFCs <- hcobject[["integrated_output"]][["GFC_all_layers"]]
  network <- hcobject[["integrated_output"]][["cluster_calc"]][["network_col_by_module"]]

  GFCs_wog <- GFCs[, base::colnames(GFCs)[!base::colnames(GFCs) == "Gene"]]
  if(base::is.vector(GFCs_wog)){
    GFCs_wog <- base::data.frame(GFCs_wog)
  }
  colors <- base::apply(GFCs_wog,2, function(x){
    y <- ((base::round(x, digits = 1) +2) *10)+1
    grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(base::length(base::seq(-2, 2, by = .1)))[y] 
  }) %>% base::as.data.frame() 
  colors$name <- GFCs$Gene
  colors[] <- base::lapply(colors, as.character)
  
  gene_to_cluster <- base::do.call(rbind, base::apply(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], 1, function(x){
    tmp <- x["gene_n"] %>%
      base::strsplit(., split = ",")%>%
      base::unlist(.)
    base::data.frame(name = tmp, color = base::rep(x["color"], base::length(tmp)))
  }))
  
  gene_to_cluster$name <- base::as.character(gene_to_cluster$name)
  gene_to_cluster$color <- base::as.character(gene_to_cluster$color)
  if(base::length(igraph::V(network)$name[!igraph::V(network)$name %in% gene_to_cluster$name]) > 0){
    gene_to_cluster <- base::rbind(gene_to_cluster, base::data.frame(name = igraph::V(network)$name[!igraph::V(network)$name %in% gene_to_cluster$name],
                                                         color = "white"))
    gene_to_cluster$name <- base::as.character(gene_to_cluster$name)
    gene_to_cluster$color <- base::as.character(gene_to_cluster$color)
  }
  colors <- colors[base::match(igraph::V(network)$name, colors$name),]
  
  if(base::any((gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name"))%in% 
         igraph::get.vertex.attribute(network)$name)){
    remove_v <- gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name")
    remove_v <- remove_v[remove_v %in% igraph::V(network)$name]
    network <- igraph::delete.vertices(network, remove_v)
  }
  
  colors <- colors[!colors$name %in% (gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name")), ]
  l <- hcobject[["integrated_output"]][["cluster_calc"]][["layout"]][base::match(igraph::V(network)$name, base::rownames(hcobject[["integrated_output"]][["cluster_calc"]][["layout"]])), ]
  for(x in base::colnames(colors)){
    if(! x == "name"){
      igraph::V(network)$color <- colors %>% dplyr::pull(., var = x)
      Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Network_GFC_", x, ".pdf"), 
                      width = 15, height = 15)
      igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3,
                          layout = l,
                          main = x)
      dev.off()
      igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3,
                          layout = l, 
                          main = x)
    }
    
  }
}