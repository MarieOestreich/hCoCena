#' Cluster Calculation
#' 
#' The function offers several community detection algorithms to identify dense regions in the co-expression network.
#'  These dense regions represent collections of highly co-expressed genes likely to form a functional group.
#' @param cluster_algo The clustering algorithm to be used. The choice is between "cluster_leiden" (default), "cluster_louvain", "cluster_label_prop", "cluster_fast_greedy",
#'  "cluster_infomap", "cluster_walktrap" and "auto" (in which case all are tested and the one with the highest modularity is chosen).
#' @param no_of_iterations Some of the algorithms are iterative (e.g. Leiden Algorithm). Set here, how many iterations should be performed. 
#'  For information on which other algorithms are iterative, please refer to their documentation in the igraph or leidenbase package. Default is 2.
#'  @param max_cluster_count_per_gene The maximum number of different clusters a gene is allowed to be associated with during the different iterations before it is marked as indecisive and removed.
#'  Default is 1.
#' @param resolution The cluster resolution if the cluster algorithm is set to "cluster_leiden". Default is 0.1. Higher values result in more clusters and vice versa.
#' @param partition_type Name of the partition type. Select from 'CPMVertexPartition', 'ModularityVertexPartition', 'RBConfigurationVertexPartition' and 'RBERVertexPartition' Default is 'ModularityVertexPartition.

#' @export 

cluster_calculation <- function(cluster_algo = "cluster_leiden",
                                no_of_iterations = 2,
                                resolution = 0.1,
                                partition_type = "ModularityVertexPartition",
                                max_cluster_count_per_gene = 1,
                                return_result = F) {
  
  hcobject[["global_settings"]][["chosen_clustering_algo"]] <<- cluster_algo
  
  if(cluster_algo == "cluster_leiden"){
    hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]] <<- leiden_clustering(g = hcobject[["integrated_output"]][["merged_net"]], num_it = no_of_iterations, resolution = resolution, partition_type = partition_type)
  }else{
    
    # count the number of graph components present in the network:
    comps <- igraph::count_components(hcobject[["integrated_output"]][["merged_net"]])
    
    # Cluster algorithms they are available:
    cluster_algo_list <- c("cluster_label_prop",
                           "cluster_fast_greedy",
                           "cluster_louvain",
                           "cluster_infomap",
                           "cluster_walktrap",
                           "cluster_leiden")
    
    
    if(cluster_algo == "auto"){
      algos_to_use <- cluster_algo_list
    }else{
      algos_to_use <- cluster_algo
    }
    
    set.seed(168575)
    
    if(cluster_algo == "auto"){
      print(algos_to_use)
      df_modularity_score <- base::do.call("rbind", base::lapply(algos_to_use,
                                                                 cluster_calculation_internal,
                                                                 graph_obj = hcobject[["integrated_output"]][["merged_net"]],
                                                                 resolution = resolution,
                                                                 partition_type = partition_type,
                                                                 case = "test"))
      
      
      cluster_algo_used <- df_modularity_score %>%
        dplyr::filter(modularity_score == base::max(modularity_score)) %>%
        dplyr::select(cluster_algorithm) %>%
        base::as.character()
      
      print(base::paste(cluster_algo_used, " will be used based on the highest modularity score."))
      
    }else{
      
      cluster_algo_used <- cluster_algo
      print(paste(cluster_algo_used, " will be used based on your input."))
      
    }
    
    
    igraph_list <- list()
    igraph_list[[1]] <- hcobject[["integrated_output"]][["merged_net"]]
    
    # apply best clustering algorithm:
    gene_which_cluster <- base::do.call("cbind", base::lapply(1:no_of_iterations,
                                                              cluster_calculation_internal,
                                                              algo = cluster_algo_used,
                                                              case = "best",
                                                              resolution = resolution,
                                                              partition_type = partition_type,
                                                              graph_obj = hcobject[["integrated_output"]][["merged_net"]]))
    
    
    
    # frequency and identity of cluster assingment of genes
    if(base::ncol(gene_which_cluster) > 1) {
      gene_cluster_ident <- base::apply(gene_which_cluster, 1, function(x){
        if(base::length(base::unique(x)) > max_cluster_count_per_gene) {   
          0
        }else{
          base::names(base::which(base::table(x) == base::max(base::table(x))))[1]
        }
      })
    } else {
      gene_cluster_ident <- gene_which_cluster[,1]
    }
    
    
    
    white_genes_clustercounts <- base::as.integer(base::grep(gene_cluster_ident, pattern = "\\b0\\b") %>%
                                                    base::length() %>% base::as.character())
    
    
    print(base::paste(white_genes_clustercounts, "genes were assigned to more than", max_cluster_count_per_gene,
                      "cluster(s). These genes are assigned to Cluster 0 (white) and will be left out of the network and further analyses."))
    
    
    
    cluster_Data <- base::data.frame(genes = igraph::vertex_attr(hcobject[["integrated_output"]][["merged_net"]], "name"),
                                     clusters = base::paste0("Cluster ", gene_cluster_ident),
                                     stringsAsFactors = FALSE)
    
    #summarize the data
    #produces a table where col are cluster name, number of components,
    #names of genes in cluster
    
    dfk <- cluster_Data %>%
      dplyr::count(clusters,genes) %>%
      dplyr::group_by(clusters) %>%
      dplyr::summarise(gene_no = base::sum(n), gene_n = base::paste0(genes, collapse = ",")) %>%
      dplyr::mutate(cluster_included = base::ifelse(gene_no >= hcobject[["global_settings"]][["min_nodes_number_for_cluster"]], "yes", "no"), color = "white")
    
    
    ##ggplot
    
    color.cluster <- get_cluster_colours()
    
    plot_clusters <- ggplot2::ggplot(data = dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", ], ggplot2::aes(x = clusters)) +
      ggplot2::geom_bar(ggplot2::aes(fill = clusters)) +
      ggplot2::scale_fill_manual(values = color.cluster)
    
    plot_clust <- ggplot2::ggplot_build(plot_clusters)
    
    dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", "color" ] <- plot_clust$data[[1]]["fill"]
    
    
    
    
    white_genes_clustersize <- base::as.integer(dfk %>% dplyr::filter(cluster_included=="no") %>%
                                                  dplyr::summarise(n = base::sum(gene_no)) %>% 
                                                  purrr::map(1))
    
    
    print(base::paste0(white_genes_clustersize, " genes were assigned to clusters with a smaller size than the defined minimal cluster size of ",
                       hcobject[["global_settings"]][["min_nodes_number_for_cluster"]], " genes per cluster. These genes will also be assigned to Cluster 0 (white) and left out of the network and further analyses."))
    
    
    
    
    # for each cluster produces a row of means per condition (info data) for all genes within the cluster
    
    gfc_dat <- hcobject[["integrated_output"]][["GFC_all_layers"]]
    
    dfk_allinfo <- base::do.call("rbind", base::lapply(1:base::nrow(dfk), gfc_mean_clustergene, 
                                                       cluster_df = dfk,
                                                       gfc_dat = hcobject[["integrated_output"]][["GFC_all_layers"]]))
    
    dfk_allinfo$vertexsize <- base::ifelse(dfk_allinfo$cluster_included == "yes", 3, 1)
    
    if(return_result == T){
      return(dfk_allinfo)
    }else{
      hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]] <<- dfk_allinfo
    }
  }
  
}
