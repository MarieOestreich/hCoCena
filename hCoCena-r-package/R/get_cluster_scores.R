#' Get Cluster Scores
#' 
#' For every gene, the ratio of its edges to genes in the same cluster to its total number of edges is determined. 
#' 	The corresponding values are returned as a data frame and a box plot is generated showing the scores for each of the clusters.
#' @param save A Boolean. Whether or not to save the plot to PDF (default is TRUE).
#' @export

get_cluster_scores <- function(save = T){

  gtc <- GeneToCluster()
  # remove white cluster since not of interest:
  gtc <- dplyr::filter(gtc, !color == "white")

  e <- hcobject[["integrated_output"]][["combined_edgelist"]]

  gtc$score <- base::lapply(gtc$gene, function(x){
  	# filter all edges that contain current gene:
    tmp <- e[e$V1 == x | e$V2 == x,]
    # get number of edges that end in current gene:
    nedges <- base::nrow(tmp)

    # get cluster colour of gene:
    c <- dplyr::filter(gtc, gene == x) %>% dplyr::pull(., "color")
    # get neighbours of current gene:
    neighbours <- base::unique(base::c(tmp$V1, tmp$V2))
    neighbours <- neighbours[!neighbours == x]
    # get neighbouring genes that are in the same cluster:
    neighbour_colours <- dplyr::filter(gtc, gene %in% neighbours & color == c)
    # get number of edges the current gene has to other genes in its own cluster:
    nedges_in_cluster <- base::nrow(neighbour_colours)
    # ration of intra-clsuter edges to total edges of current gene:
    return(nedges_in_cluster/nedges)
  }) %>% base::unlist()

  
  gtc$label <- base::lapply(gtc$gene, function(x){
    c <- dplyr::filter(gtc, gene == x) %>% dplyr::pull(., color)
    tmp <- dplyr::filter(gtc, color == c)
    return(base::paste0(c, " [", base::nrow(tmp), "]"))
  }) %>% base::unlist()


  p <- ggplot2::ggplot(gtc, ggplot2::aes(x=label, y= score, color=color, fill = color)) +
    ggplot2::geom_boxplot()+
    ggplot2::scale_color_manual(values = base::sort(base::unique(gtc$color))) +
    ggplot2::scale_fill_manual(values = grDevices::adjustcolor(base::sort(base::unique(gtc$color)), alpha.f = 0.5)) +
    ggplot2::xlab("cluster")+
    ggplot2::theme_bw()+ggplot2::coord_flip() + 
    ggplot2::ggtitle("Cluster scores")

  graphics::plot(p)

  if(save){
    Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Cluster_scores.pdf"), 
                    width = 8)
    
    graphics::plot(p)
    
    grDevices::dev.off()
  }

  gtc$label <- NULL
  hcobject[["satellite_outputs"]][["cluster_scores"]] <<- list(scores_per_gene = gtc, plot = p)
  
}