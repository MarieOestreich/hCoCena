#' Plot Integrated Network
#' 
#' Creates an igraph of the integrated network.
#' @param layout NULL or a precalculated layout. If NULL (default) the layout is calculated.
#' @param gene_labels A vector of character strings containing the names of genes to be labelled. Labels will be placed with an offset, which can be adjusted using the label_offset parameter.
#' @param save A Boolean stating whether of not to save the network to a PDF. Default is TRUE.
#' @param label_offset The offset of provided gene labels to the network center. Default is 50, depends on network size.
#' @export

plot_integrated_network <- function(layout = NULL, 
									gene_labels = NULL, 
									save = T, 
									label_offset = 50){

  network <- hcobject[["integrated_output"]][["merged_net"]]

  if(hcobject[["global_settings"]][["layout_algorithm"]] == "cytoscape"){
    if(is.null(layout)){
      message("Please calculate the Cytoscape layout in the satellite markdown. If you have already done that, set the 'layout' parameter to hcobject[['integrated_output']][['cluster_calc']][['layout']].")
      return()
    }
  }
  gene_to_cluster <- GeneToCluster()
  colnames(gene_to_cluster) <- c("name", "color")
 
  gene_to_cluster$name <- base::as.character(gene_to_cluster$name)
  gene_to_cluster$color <- base::as.character(gene_to_cluster$color)
  
  if(base::length(igraph::V(network)$name[!igraph::V(network)$name %in% gene_to_cluster$name]) > 0){
    gene_to_cluster <- base::rbind(gene_to_cluster, base::data.frame(name = igraph::V(network)$name[!igraph::V(network)$name %in% gene_to_cluster$name],
                                                         color = "white"))
    gene_to_cluster$name <- base::as.character(gene_to_cluster$name)
    gene_to_cluster$color <- base::as.character(gene_to_cluster$color)
  }

  gene_to_cluster <- gene_to_cluster[base::match(igraph::V(network)$name, gene_to_cluster$name),]

  igraph::V(network)$color <- gene_to_cluster$color
  if("white" %in% base::unique(gene_to_cluster$color)){
    network <- igraph::delete.vertices(network, gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name"))
  }
  

  if(base::is.null(layout)){
    base::set.seed(123)
    if(hcobject[["global_settings"]][["layout_algorithm"]] == "cytoscape"){
  
      if("layout" %in% base::names(hcobject[["integrated_output"]][["cluster_calc"]])){
        l <- base::as.matrix(hcobject[["integrated_output"]][["cluster_calc"]][["layout"]])
      }else{
        stop("Please create a Cytoscape-based network layout first. To do so, refer to the satellite markdown, section 'Cytoscape'")
      }
    }else{
      l <- igraph::layout.fruchterman.reingold(network)
      base::rownames(l) <- igraph::V(network)$name
    }
    
  }else{
    l <- layout %>% base::as.data.frame()
    l <-l[base::match(igraph::get.vertex.attribute(network)$name, base::rownames(l)),]
    l[] <- base::lapply(l, as.numeric)
    l <- base::as.matrix(l)
  }
  
  if(!base::is.null(gene_labels)){
    
    new_genes_df <- base::data.frame(name = base::paste0("label_", base::c(1:base::length(gene_labels))), label = gene_labels %>% base::as.character())
    mean_coord_x <- stats::median(l[,1])
    mean_coord_y <- stats::median(l[,2])
    new_indeces <- base::match(new_genes_df$label, igraph::get.vertex.attribute(network)$name)
    new_genes_df$coords_x <- l[new_indeces,1]
    new_genes_df$coords_y <- l[new_indeces,2]
    left_up <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y >= mean_coord_y)%>%
      dplyr::pull(., label)
    left_down <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y < mean_coord_y)%>%
      dplyr::pull(., label)
    right_up <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y >= mean_coord_y)%>%
      dplyr::pull(., label)
    right_down <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y < mean_coord_y)%>%
      dplyr::pull(., label)
    
    new_position_l <- base::matrix(base::cbind(base::rep(base::ceiling(base::min(l[,1]))-label_offset, base::length(left_up)+ base::length(left_down)), 
                                 base::seq(from=base::ceiling(base::max(l[,2])), to = base::ceiling(base::min(l[,2])), 
                                     length.out = base::length(left_up)+ base::length(left_down))), ncol = 2)
    new_position_r <- base::matrix(base::cbind(base::rep(base::ceiling(base::max(l[,1]))+label_offset, base::length(right_up)+ base::length(right_down)), 
                                   base::seq(from=base::ceiling(base::max(l[,2])), to = base::ceiling(base::min(l[,2])), 
                                       length.out = base::length(right_up)+ base::length(right_down))), ncol = 2)
    new_position <- base::rbind(new_position_l, new_position_r)
    
    base::colnames(new_position) <- base::colnames(l)
    l2 <- base::rbind(l, new_position) %>%
    	base::as.matrix()
    new_genes_df_l <- new_genes_df[new_genes_df$label %in% base::c(left_up, left_down),]
    new_genes_df_l <- new_genes_df_l[base::order(new_genes_df_l$coords_y, decreasing = T),]
    new_genes_df_r <- new_genes_df[new_genes_df$label %in% base::c(right_up, right_down),]
    new_genes_df_r <- new_genes_df_r[base::order(new_genes_df_r$coords_y, decreasing = T),]
    new_genes_df <- base::rbind(new_genes_df_l, new_genes_df_r)
    
    
    network2 <- igraph::add.vertices(network, nv = base::length(new_genes_df$name), 
                            attr = list(name = new_genes_df$name))
    new_labels <- base::lapply(igraph::get.vertex.attribute(network2)$name, function(x){
      if(x %in% igraph::get.vertex.attribute(network)$name){
        NA
      }else{
        dplyr::filter(new_genes_df, name == x)%>%
          dplyr::pull(., "label")
      }
    })%>% base::unlist()
    
    igraph::V(network2)$label <- new_labels
    
    new_edges <- base::matrix(base::c(base::match(new_genes_df$label, igraph::get.vertex.attribute(network2)$name), 
                        base::match(new_genes_df$name, igraph::get.vertex.attribute(network2)$name)), 
                        ncol = 2, byrow = F)
    new_edges <- base::as.vector(base::t(new_edges))
      
    network2 <- igraph::add.edges(network2, new_edges)
    
    new_edge_color <- base::apply(igraph::get.edgelist(network2),1, function(x){
      if(x[2] %in% new_genes_df$name){
        igraph::get.vertex.attribute(network2, name = "color", index = x[1])
      }else{
        "lightgrey"
      }
    })
    
    new_label_color <- base::lapply(igraph::get.vertex.attribute(network2)$name, function(x){
      if(x %in% new_genes_df$name){
        tmp_gene_n <- dplyr::filter(new_genes_df, name == x)%>%
          dplyr::pull(., label)
        igraph::get.vertex.attribute(network2, name = "color", index = tmp_gene_n)
      }else{
        NA
      }
    })%>% base::unlist()
    
    
    new_frame_color <- base::lapply(igraph::get.vertex.attribute(network2)$name, function(x){
      if(x %in% new_genes_df$name){
        "white"
      }else{
        "black"
      }
    })%>% base::unlist()
    
    if(save == T){
      Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Network_col_by_module.pdf"), 
                      width = 15, height = 15)
      
      igraph::plot.igraph(network2, vertex.size = 3, vertex.label = new_labels, vertex.label.cex = 0.75,
                          layout = l2, main = "Co-expression network coloured by module",
                          edge.color = new_edge_color, vertex.label.color = new_label_color,
                          vertex.frame.color = new_frame_color)
      
      grDevices::dev.off()
    }
    
    
    igraph::plot.igraph(network2, vertex.size = 3, vertex.label = new_labels, vertex.label.cex = 0.75,
                        layout = l2, main = "Co-expression network coloured by module",
                        edge.color = new_edge_color, vertex.label.color = new_label_color,
                        vertex.frame.color = new_frame_color)
    
    hcobject[["integrated_output"]][["cluster_calc"]][["labelled_network"]] <<- network2
    hcobject[["integrated_output"]][["cluster_calc"]][["network_col_by_module"]] <<- network
  }else{
    
    if(save == T){
      Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Network_modules.pdf"), 
                      width = 15, height = 15)
      
      igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3, 
                          layout = l, main = "Co-expression network coloured by module")
      
      grDevices::dev.off()
    }
    
    
    igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3, 
                        layout = l, main = "Co-expression network coloured by module")
    
    hcobject[["integrated_output"]][["cluster_calc"]][["labelled_network"]] <<- list()
    hcobject[["integrated_output"]][["cluster_calc"]][["network_col_by_module"]] <<- network
  }
  
  
  if(base::is.null(layout)){
    hcobject[["integrated_output"]][["cluster_calc"]][["layout"]] <<- l
  }
  
}