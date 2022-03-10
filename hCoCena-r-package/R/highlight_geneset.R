#' Highlight Geneset
#' 
#' Receives and highlights a list of genes in the integrated network.
#' @param gene_set A vector of gene symbols (as strings) to be highlighted in the network.
#' @param col A string specifying the color with which the nodes of the genes are framed to highlight them. Default is "black".
#' @export

highlight_geneset <- function(gene_set, col = "black"){
  
  gtc <- GeneToCluster() %>% dplyr::filter(., !color == "white")
  g <- hcobject[["integrated_output"]][["merged_net"]]
  g <- igraph::delete_vertices(g, igraph::V(g)$name[!igraph::V(g)$name %in% gtc$gene])
  
  base::rownames(gtc) <- gtc$gene
  gtc <- gtc[igraph::V(g)$name,]
  
  igraph::V(g)$color <- base::lapply(gtc$gene, function(x){
    if(x %in% gene_set){
      dplyr::filter(gtc, gene == x) %>% dplyr::pull(., "color")
    }else{
      "white"
    }
  }) %>% base::unlist()
  
  igraph::V(g)$size<- base::lapply(gtc$gene, function(x){

    if(x %in% gene_set){
      5
    }else{
      3
    }
  }) %>% base::unlist()
  
  
  igraph::V(g)$frame.color <- base::lapply(gtc$gene, function(x){
    #c <- dplyr::filter(gtc, gene == x) %>% dplyr::pull(., "color")
    if(x %in% gene_set){
      col
    }else{
      "lightgrey"
    }
  }) %>% base::unlist()
  
  l <- hcobject[["integrated_output"]][["cluster_calc"]][["layout"]]
  l <- l[base::rownames(l) %in% igraph::V(g)$name,]
  l <- l[igraph::V(g)$name, ]
  igraph::plot.igraph(g, vertex.label = NA, layout = l)
  
}