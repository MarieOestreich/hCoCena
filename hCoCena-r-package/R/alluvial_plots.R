#' Plot Alluvial Plots For Different Clustering Algorithms
#' 
#' Alluvial plots will be generated demonstrating how the genes change clusters if the current clustering algorithm (always shown on the left) was changed to any of the other clustering options (always shown on the right).
#' 	The produced plot is interactive, hover over it to get more details.
#' @export

algo_alluvial <- function(){

  network <- hcobject[["integrated_output"]][["merged_net"]]

  cluster_algo_list <- base::c("cluster_louvain",
                       "cluster_fast_greedy",
                       "cluster_infomap",
                       "cluster_walktrap",
                       "cluster_label_prop",
                       "cluster_leiden")


  color.cluster <- get_cluster_colours()


  current_algo <- NULL

  for(c in base::unique(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]$color)){
    genes <- dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], color == c) %>%
      dplyr::pull(., "gene_n") %>%
      base::strsplit(., split = ",") %>%
      base::unlist(.)
    current_algo <- base::rbind(current_algo, base::data.frame(gene = genes, cluster = base::rep(c, base::length(genes))))
  }

  current_algo$cluster <- base::as.factor(current_algo$cluster)

  output <- list()
  output[[hcobject[["global_settings"]][["chosen_clustering_algo"]]]] <- current_algo
  for(a in cluster_algo_list[!cluster_algo_list == hcobject[["global_settings"]][["chosen_clustering_algo"]]]){

    base::set.seed(1)
    base::set.seed(.Random.seed[1])

    # add exception for leiden:
    if(a == "cluster_leiden"){
      partition <- leidenAlg::leiden.community(graph = network)
      partition_df <- base::data.frame(gene = partition$names, cluster = base::as.numeric(partition$membership))
      partition_df$cluster <- partition_df$cluster + 1
    }else{
      cfg = base::get(a)(network)
      partition_df <- base::data.frame(gene = igraph::get.vertex.attribute(network)$name, cluster = cfg$membership)
      
    }

    partition_df$cluster <- base::lapply(partition_df$cluster, function(x){
      color.cluster[x]
    })%>% base::unlist()
    
    
    
    partition_df_table <- base::table(partition_df$cluster) %>% 
    	base::as.data.frame()
    partition_df <- partition_df[partition_df$cluster %in% (dplyr::filter(partition_df_table, Freq >= hcobject[["global_settings"]][["min_nodes_number_for_cluster"]]) %>%
                 dplyr::pull(., "Var1")),]
    partition_df$cluster <- base::as.factor(partition_df$cluster)


    a2a <- base::merge(current_algo, partition_df , by = "gene", all = T)
    base::colnames(a2a) <- base::c("gene", "cluster1", "cluster2")
    a2a$cluster1 <- base::as.character(a2a$cluster1)
    a2a$cluster2 <- base::as.character(a2a$cluster2)
    a2a[base::is.na(a2a)] <- "white"
    a2a$Freq <- 1
    a2a$merged <- base::paste0(a2a$cluster1, a2a$cluster2)



    links <- a2a[, 2:base::ncol(a2a)] %>% 
    	base::unique()
    links$Freq <- base::lapply(links$merged, function(x){
      base::length(a2a$merged[a2a$merged == x])
    }) %>% 
    	base::unlist()
    links$merged <- NULL
    links$cluster1_name <- base::paste0(links$cluster1, "_", base::strsplit(hcobject[["global_settings"]][["chosen_clustering_algo"]], split = "_")[[1]][2])
    links$cluster2_name <- base::paste0(links$cluster2, base::strsplit(a, split = "cluster")[[1]][2])

    nodes <- base::data.frame(name = base::c(base::as.character(links$cluster1_name) %>% 
    		base::unique(), base::as.character(links$cluster2_name) %>% 
    		base::unique()))

    links$IDsource <- base::match(links$cluster1_name, nodes$name)-1
    links$IDtarget <- base::match(links$cluster2_name, nodes$name)-1
    base::colnames(links) <- base::c("source_col", "target_col", "value", "source", "target","IDsource", "IDtarget")
    node_col <- base::data.frame(name=base::c(base::as.character(links$source_col)%>% base::unique(), base::as.character(links$target_col)%>% base::unique()))
    my_color <- base::paste0('d3.scaleOrdinal() .domain([', base::paste0(base::paste0('"', base::as.character(nodes$name), collapse = '", '),'"', collapse = '"'),
                       ']) .range([', base::paste0(base::paste0('"', base::as.character(node_col$name), collapse = '", '), '"',collapse = '"'), '])')
    p <- networkD3::sankeyNetwork(Links = links, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name",
                       sinksRight=FALSE, colourScale=my_color,
                       nodeWidth=40, fontSize=13, nodePadding=10)
  
    print(p)
    output[[a]] <- partition_df
  }

  hcobject[["integrated_output"]][["alluvials"]] <<- output
  
}
