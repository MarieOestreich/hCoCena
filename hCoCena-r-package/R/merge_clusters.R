#' Merge Clusters
#' 
#' Merge similar clusters in module heatmap.
#' @param k An integer. The number of desired clusters created by cutting the clustering tree.
#' @param save Boolean. If FALSE (default) cutting is showcased but not saved. Set to TRUE once you have found the right clsutering and overwrite the old cluster (cannot be reversed).
#' @param method Method used for clustering. Default is "complete".
#' @export

merge_clusters <- function(k = 1, save = F, method = "complete"){
  m <- hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]]@ht_list[[1]]@matrix
  p <- pheatmap::pheatmap(mat = m,
                          color = base::rev(RColorBrewer::brewer.pal(11, "RdBu")),
                          scale = "none",
                          cluster_rows =  T,
                          cluster_cols =  T,
                          fontsize = 8,
                          show_rownames = F, 
                          show_colnames = T, 
                          clustering_distance_cols = "euclidean", 
                          clustering_method = method, 
                          cutree_rows = k)
  if(save){
    cut <- stats::cutree(p$tree_row, k=k)
    new_group <- base::data.frame(old_cluster = base::names(cut), group = cut)
    new_group$new_cluster <- base::lapply(new_group$group, function(x){
      get_cluster_colours()[x]
    }) %>% base::unlist()
    
    gtc <- GeneToCluster()
    base::colnames(gtc) <- base::c("gene", "old_color")
    gtc$color <- base::lapply(gtc$old_color, function(x){
      nc <- dplyr::filter(new_group, as.character(old_cluster) == x) %>% dplyr::pull(., "new_cluster")
      
      return(nc)
    }) %>% base::unlist()
    
    new_cluster_info <- NULL
    
    for(c in base::unique(gtc$color)){
      tmp <- gtc[gtc$color == c, ]
      gene_no <- base::nrow(tmp)
      gene_n <- base::paste0(tmp$gene, collapse = ",")
      if(c == "white"){
        cluster_included <- "no"
        vertexsize <- 1
      }else{
        cluster_included <- "yes"
        vertexsize <- 3
      }
      
      color <- c
      conditions <- base::paste0(base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])[1:(base::ncol(hcobject[["integrated_output"]][["GFC_all_layers"]])-1)], collapse = "#")
      gfc_means = hcobject[["integrated_output"]][["GFC_all_layers"]][hcobject[["integrated_output"]][["GFC_all_layers"]][["Gene"]] %in% tmp$gene,] %>%
        dplyr::select(-Gene) %>%
        base::colMeans()
      
      grp_means = base::paste0(base::round(gfc_means,3) , collapse = ",")
      
      new_cluster_info <- rbind(new_cluster_info,
                                data.frame(clusters = c,
                                           gene_no = gene_no,
                                           gene_n = gene_n,
                                           cluster_included = cluster_included,
                                           color = c,
                                           conditions = conditions,
                                           grp_means = grp_means,
                                           vertexsize = vertexsize,
                                           stringsAsFactors = F))
      
    }
    hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]] <<- new_cluster_info
    plot_cluster_heatmap()
    
  }else{
    p
    message("This is only a preview. To save these new clusters, set the save parameter to TRUE.")
  }
  
}