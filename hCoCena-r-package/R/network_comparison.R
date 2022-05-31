#' Compare clusters of 2 networks
#' 
#' This function calculates and visualizes the Jaccard-Index of all pairs of clusters from two networks.
#' This allows the comparison of the two networks with respect to the clusters they form and how those clusters relate to each other.
#' @param gtc1_path File path to a 'gtc'-file (gene-to-clsuter) of the first network: Essentially just a file with two columns, the first containing gene symbols and the second 
#' giving the cluster each gene belongs to. Such a file can be generated during an hCoCena analysis using the function 'hcocena::export_clusters()', but as long as
#' the described structure is preserved it can also be generated manually elsewhere.
#' @param gtc2_path See 'gtc1_path', only for network 2.
#' @param sep The separator of the 'gtc'-file. Default is tab-separated.
#' @param header A Boolean. Whether or not the file has headers (column names).
#' @param cellsize The size of the cells/tiles in the plotted heatmap. Default is 18, may be adjusted for aestetics reasons.
#' @return The heatmap object for replotting/re-sizing etc. and the result matrix. Output can be found under hcobject$satellite_outputs$network_comparison_1
#' @export

network_comparison_1 <- function(gtc1_path, gtc2_path, sep = "\t", header = TRUE, cellsize = 18){
  
  gtc1 <- readr::read_delim(file = gtc1_path , delim = sep, col_names = header)
  base::colnames(gtc1) <- c('gene', 'color')
  gtc2 <- readr::read_delim(file = gtc2_path , delim = sep, col_names = header)
  base::colnames(gtc2) <- c('gene', 'color')
  
  out <- base::list()
  
  for(c1 in unique(gtc1$color)){
    cdf <- NULL

    # get gene set of this cluster:
    set1 <- dplyr::filter(gtc1, color == c1) %>% dplyr::pull(., 'gene')
    
    # for plot indicate which network this cluster came from:
    cname1 <- base::paste0(c1, ' network 1 [',length(set1), ']')
    
    for(c2 in unique(gtc2$color)){
      
      # get gene set of this cluster:
      set2 <- dplyr::filter(gtc2, color == c2) %>% dplyr::pull(., 'gene')
      # for plot indicate which network this cluster came from:
      cname2 <- base::paste0(c2, ' network 2 [',length(set2), ']')
      # calculate Jaccard Index (JI):
      JI <- calc_jaccard(set1, set2)
      
      rdf <- base::data.frame(V1 = JI)
      base::colnames(rdf) <- c(cname1)
      base::rownames(rdf) <- c(cname2)
      cdf <- base::rbind(cdf, rdf)
    }
    
    out[[c1]] <- cdf
    
  }
  out <- rlist::list.cbind(out)
  p <- pheatmap::pheatmap(as.matrix(out), 
                     cluster_rows = F, 
                     cluster_cols = F, 
                     cellheight = cellsize, 
                     cellwidth = cellsize,
                     color = RColorBrewer::brewer.pal(name='Blues', n=9),
                     display_numbers = T, 
                     fontsize_number = 6, 
                     number_color = 'orange', main = 'Jaccard Index of Cluster Pairs')
  # save to PDF
  Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/network_comparison_1.pdf"), 
                  width = 10, 
                  height = 7)
  p
  grDevices::dev.off()
  
  hcobject[["satellite_outputs"]][["network_comparison_1"]] <<-list(heatmap=p, matrix=as.matrix(out))
}


#' Calculate Jaccard Index of two sets.
#' 
#' Only for internal use.
#' @noRd

calc_jaccard <- function(set1, set2){
  # get size of the intersection of the 2 sets:
  intersection <- base::intersect(set1, set2) %>% length()
  # get size of the union of the 2 sets:
  union <- base::union(set1, set2) %>% length()
  # return the jaccard index of the 2 sets (|intersection|/|union|):
  return(intersection/union)
}




#' Compare degree distribution of gene set
#' 
#' This function accepts two networks and a set of genes. It then calculates the degree-distribution of each gene in both networks and also the Jaccard-Index
#' of each gene's neighbourhoods in the two networks. The results are visualized in a 2D dot plot.
#' This allows the comparison of the two networks with respect to the connectivity of selected genes.
#' @param net1 An igraph object of the first network. hCoCena provides its integrated network as igraph objects. You can find it here: hcobject$integrated_output$merged_net.
#' If you want to create a network manually with external data, please refer to the igraph documentation: https://igraph.org/r/.
#' @param net2 See 'net1', only for network 2.
#' @param as In the future we plan to provide either an igraph object or an edge list. For now, only the igraph option is available.
#' @param gene_vec A vector of gene symbols you wish to investigate. Gene symbols must be provided as strings, e.g., c( "YME1L1", "SLC2A5", "SIAH2", "GPI", "IL10RB").
#' @return The ggplot object for replotting/re-sizing/modification etc. and the data used to plot the ggplot. Output can be found under hcobject$satellite_outputs$network_comparison_2
#' @export

network_comparison_2 <- function(net1, net2, as = 'igraph', gene_vec){
  nodes1 <- igraph::V(net1)$name
  nodes1 <- nodes1[nodes1 %in% gene_vec]
  nodes2 <- igraph::V(net2)$name
  nodes2 <- nodes2[nodes2 %in% gene_vec]
  
  intersection <- base::intersect(nodes1, nodes2)
  if(length(intersection) == 0){
    base::message('None of the given genes are present in both networks. Aborting.')
    return(NULL)
  }else{
    base::message(length(intersection), ' of the given genes are present in both networks.')
  }
  
  out <- NULL
  for(gene in intersection){
    # neighbours in net1:
    neighbours1 <- igraph::neighbors(net1, gene, mode = c("all"))$name
    # neighbours in net2:
    neighbours2 <- igraph::neighbors(net2, gene, mode = c("all"))$name
    JI_common = (base::intersect(neighbours1, neighbours2) %>% length())/(base::union(neighbours1, neighbours2) %>% length())
    rdf <- base::data.frame(gene = gene, 
                            neighbours_in_1 = log(length(neighbours1), base = 2), 
                            neighbours_in_2 = log(length(neighbours2), base = 2), 
                            JI_neighbours = JI_common)
    out <- base::rbind(out, rdf)
  }
  
  g <- ggplot2::ggplot(out, ggplot2::aes(x=neighbours_in_1, y=neighbours_in_2,label=gene))+
    ggplot2::geom_point(ggplot2::aes(x=neighbours_in_1, y=neighbours_in_2, size=JI_neighbours, color=JI_neighbours))+
    ggplot2::theme_bw()+
    ggplot2::geom_text(hjust=0.5, vjust=-1, size=3)+
    ggplot2::ylim(c(min(c(out$neighbours_in_1, out$neighbours_in_2)), max(c(out$neighbours_in_1, out$neighbours_in_2))))+
    ggplot2::xlim(c(min(c(out$neighbours_in_1, out$neighbours_in_2)), max(c(out$neighbours_in_1, out$neighbours_in_2))))+
    ggplot2::ylab("node degree in network 2 (log2)")+
    ggplot2::xlab("node degree in network 1 (log2)")+
    ggplot2::geom_abline()+
    ggplot2::guides(color = ggplot2::guide_colorbar(order = 1),
                    size = ggplot2::guide_legend(order = 2))+
    ggplot2::labs(color = "Jaccard-Index of Neighbours", size=ggplot2::element_blank())
  graphics::plot(g)
  hcobject[["satellite_outputs"]][["network_comparison_2"]] <<- list(plot=g, data=out)
}




  


                                                    