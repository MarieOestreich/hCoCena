#' Highlight Geneset
#' 
#' Receives and highlights a list of genes in the integrated network.
#' @param gene_set A vector of gene symbols (as strings) to be highlighted in the network.
#' @param col A string specifying the color with which the nodes of the genes are framed to highlight them. Default is "black".
#' @param label_offset A numeric specifying the distance between the labeled genes and the network. Default is 3.
#' @param save A Boolean. Whether or not a labelled hub-network per cluster and the expression heatmap are to be save to PDF. Default is FALSE.
#' @param plot A Boolean. Wheather or not to plot the network (per cluster) with highlighted hub nodes. Default is FALSE.
#' @export

highlight_geneset <- function(gene_set, name = NULL, col = "black", label_offset = 3, plot = T, save = T){
  
  label_offset = label_offset
  g <- hcobject[["integrated_output"]][["merged_net"]]
  gtc <- GeneToCluster()
  colnames(gtc) <- c("name", "color")
  gtc$name <- base::as.character(gtc$name)
  gtc$color <- base::as.character(gtc$color)
  
  if(base::length(igraph::V(g)$name[!igraph::V(g)$name %in% gtc$name]) > 0){
    gtc <- base::rbind(gtc, base::data.frame(name = igraph::V(g)$name[!igraph::V(g)$name %in% gtc$name],
                                                                     color = "white"))
    gtc$name <- base::as.character(gtc$name)
    gtc$color <- base::as.character(gtc$color)
  }
  
  gtc <- gtc[base::match(igraph::V(g)$name, gtc$name),]
  
  igraph::V(g)$color <- gtc$color
  if("white" %in% base::unique(gtc$color)){
    g <- igraph::delete.vertices(g, gtc[gtc$color == "white",] %>% dplyr::pull(., "name"))
  }

  igraph::V(g)$color <- base::lapply(igraph::V(g)$name, function(x){
    if(x %in% gene_set){
      dplyr::filter(gtc, name == x) %>% dplyr::pull(., color)
    }else{
      NA
    }
  }) %>% base::unlist()
  
  igraph::V(g)$size <- base::lapply(igraph::V(g)$name, function(x){
    if(x %in% gene_set){
      5
    }else{
      3
    }
  }) %>% base::unlist()
  
  igraph::V(g)$frame.color <- base::lapply(igraph::V(g)$name, function(x){
    if(x %in% gene_set){
      col
    }else{
      "lightgrey"
    }
  }) %>% base::unlist()
  
  # Get original layout
  if(is.null(hcobject[["integrated_output"]][["cluster_calc"]][["layout"]])){
    if(hcobject[["global_settings"]][["layout_algorithm"]] == "cytoscape")
      stop("Please create a Cytoscape-based network layout first. To do so, refer to the satellite markdown, section 'Cytoscape'")
    message("Computing network layout using Fruchterman-Reingold algorithm...")

    l <- igraph::layout.fruchterman.reingold(g)
    base::rownames(l) <- igraph::V(g)$name
    hcobject[["integrated_output"]][["cluster_calc"]][["layout"]] <<- l
    l <- l[base::rownames(l) %in% igraph::V(g)$name,]
    l <- l[igraph::V(g)$name, ]
  } else l <- hcobject[["integrated_output"]][["cluster_calc"]][["layout"]]


  # Add genes to be highlighted (add vertices to the label and arrange plotting order)
  new_genes_df <- base::data.frame(name = base::paste0("label_", base::c(1:base::length(gene_set))), 
                                   label = gene_set %>% base::as.character())
  
  mean_coord_x <- stats::median(l[,1])
  mean_coord_y <- stats::median(l[,2])
  
  new_indeces <- base::match(new_genes_df$label, igraph::get.vertex.attribute(g)$name)
  
  new_genes_df$coords_x <- l[new_indeces, 1]
  new_genes_df$coords_y <- l[new_indeces, 2]
  
  left_up <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y >= mean_coord_y) %>% dplyr::pull(., name)
  left_down <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y < mean_coord_y) %>% dplyr::pull(., name)
  right_up <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y >= mean_coord_y) %>% dplyr::pull(., name)
  right_down <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y < mean_coord_y) %>% dplyr::pull(., name)
  
  new_position_l <- base::matrix(base::cbind(base::rep(base::ceiling(base::min(l[,1]))-label_offset, 
                                                       base::length(left_up) + base::length(left_down)), 
                                             base::seq(from = base::ceiling(base::max(l[,2])), 
                                                       to = base::ceiling(base::min(l[,2])), 
                                                       length.out = base::length(left_up)+ base::length(left_down))), 
                                 ncol = 2)
  new_position_r <- base::matrix(base::cbind(base::rep(base::ceiling(base::max(l[,1]))+label_offset, 
                                                       base::length(right_up)+ base::length(right_down)), 
                                             base::seq(from = base::ceiling(base::max(l[,2])), 
                                                       to = base::ceiling(base::min(l[,2])), 
                                                       length.out = base::length(right_up) + base::length(right_down))), 
                                 ncol = 2)
  new_position <- base::rbind(new_position_l, new_position_r)
  tmp <- c(left_up, left_down, right_up, right_down)
  base::rownames(new_position) <- tmp
  base::colnames(new_position) <- base::colnames(l)
  l2 <- base::rbind(l, new_position) %>% base::as.matrix()
  
  new_genes_df_l <- new_genes_df[new_genes_df$name %in% base::c(left_up, left_down),]
  new_genes_df_l <- new_genes_df_l[base::order(new_genes_df_l$coords_y, decreasing = T),]
  new_genes_df_r <- new_genes_df[new_genes_df$name %in% base::c(right_up, right_down),]
  new_genes_df_r <- new_genes_df_r[base::order(new_genes_df_r$coords_y, decreasing = T),]
  new_genes_df <- base::rbind(new_genes_df_l, new_genes_df_r)
  
  network2 <- igraph::add.vertices(g, nv = base::length(new_genes_df$name), attr = list(name = new_genes_df$name))

  new_label_color <- base::lapply(igraph::get.vertex.attribute(network2)$name, function(x){
    if(x %in% new_genes_df$name){
      tmp_gene_n <- dplyr::filter(new_genes_df, name == x)%>%
        dplyr::pull(., label)
      igraph::get.vertex.attribute(network2, name = "color", index = tmp_gene_n)
    }else{
      NA
    }
  })%>% base::unlist()
  
  name_to_size <- base::data.frame(name = igraph::V(network2)$name, size = igraph::V(network2)$size)
  
  vertex_size <- base::lapply(igraph::V(network2)$name, function(x){
    if(x %in% igraph::V(g)$name){
      max(l2) + abs(min(l2))
    }
    else{
      0
    }
  }) %>% base::unlist()
  
  network3 <- igraph::delete.vertices(network2, igraph::V(network2)$name[!(igraph::V(network2)$name %in% gene_set)])
  network3 <- igraph::add.vertices(network3, nv = base::length(new_genes_df$name), attr = list(name = new_genes_df$name))
  new_labels <- base::lapply(igraph::get.vertex.attribute(network3)$name, function(x){
    if(x %in% igraph::get.vertex.attribute(g)$name){
      NA
    }else{
      tmp_label <- dplyr::filter(new_genes_df, name == x) %>% dplyr::pull(., "label")
    }
  })%>% base::unlist()
  
  igraph::V(network3)$label <- new_labels
  name_to_id <- base::data.frame(name = igraph::V(network3)$name, id = 1:base::length(igraph::V(network3)$name))
  
  nti1 <- name_to_id[name_to_id$name %in% new_genes_df[,1],]
  base::rownames(nti1) <- nti1$name
  nti1 <- nti1[new_genes_df[,1],]
  
  nti2 <- name_to_id[name_to_id$name %in% new_genes_df[,2],]
  base::rownames(nti2) <- nti2$name
  nti2 <- nti2[new_genes_df[,2],]
  
  new_edges <- base::matrix(base::c(nti1$id, nti2$id), ncol = 2, byrow = F)
  
  new_edges <- base::as.vector(base::t(new_edges))
  network3 <- igraph::add.edges(network3, new_edges)
  
  igraph::V(network3)$size <- (max(l2) + abs(min(l2)))*2
  
  if(save){
    Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", name, "_network.pdf"), 
                    width = 20, height = 15)
    
    igraph::plot.igraph(network2, vertex.size = vertex_size,
                        vertex.label = NA,
                        layout = l2, 
                        edge.color = "lightgray",rescale = F, 
                        xlim = c(min(l2[,1]), max(l2[,1])), 
                        ylim = c(min(l2[,2]), max(l2[,2])))
    
    igraph::plot.igraph(network3, add = T, rescale = F, 
                        vertex.label = new_labels, 
                        edge.color = "black", 
                        vertex.label.cex = 1,
                        vertex.label.dist = 1,
                        vertex.label.color = "black",
                        layout = l2[rownames(l2) %in% igraph::V(network3)$name,], 
                        xlim = c(min(l2[,1]), max(l2[,1])), 
                        ylim = c(min(l2[,2]), max(l2[,2])))
    graphics::title(main = name, cex.main = 2)
    
    grDevices::dev.off()
  }
  if(plot){
    igraph::plot.igraph(network2, vertex.size = vertex_size,
                        vertex.label = NA,
                        layout = l2, 
                        edge.color = "lightgray",rescale = F, 
                        xlim = c(min(l2[,1]), max(l2[,1])), 
                        ylim = c(min(l2[,2]), max(l2[,2])))

    igraph::plot.igraph(network3, add = T, rescale = F, 
                        vertex.label = new_labels, 
                        edge.color = "black", 
                        vertex.label.cex = 1,
                        vertex.label.dist = 1,
                        vertex.label.color = "black",
                        layout = l2[rownames(l2) %in% igraph::V(network3)$name,], 
                        xlim = c(min(l2[,1]), max(l2[,1])), 
                        ylim = c(min(l2[,2]), max(l2[,2])))
    graphics::title(main = name, cex.main = 2)
  }
}
