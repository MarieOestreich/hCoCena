#' Colour Single Cluster
#' 
#' Replots the integrated network, highlighting a given cluster.
#' @param cluster A string, giving the colour of the cluster to be highlighted.
#' @export

colour_single_cluster <- function(cluster){

  gtc <- GeneToCluster() %>% dplyr::filter(., !color == "white")
  g <- hcobject[["integrated_output"]][["merged_net"]]
  g <- igraph::delete_vertices(g, igraph::V(g)$name[!igraph::V(g)$name %in% gtc$gene])
  
  base::rownames(gtc) <- gtc$gene
  gtc <- gtc[igraph::V(g)$name,]
  
  igraph::V(g)$color <- base::lapply(gtc$color, function(x){
    if(x == cluster){
      x
    }else{
      "white"
    }
  }) %>% base::unlist()
  igraph::V(g)$size <- base::lapply(gtc$color, function(x){
    if(x == cluster){
      5
    }else{
      3
    }
  }) %>% base::unlist()
  l <- hcobject[["integrated_output"]][["cluster_calc"]][["layout"]]
  l <- l[base::rownames(l) %in% igraph::V(g)$name,]
  l <- l[igraph::V(g)$name, ]
  igraph::plot.igraph(g, vertex.label = NA, layout = l)
  
}