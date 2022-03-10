#' Get Cluster Genes
#' 
#' Returns all genes present in a given cluster as a vector.
#' @param cluster The colour of the cluster as a string.
#' @return A vector of strings
#' @export

get_cluster_genes <- function(cluster){
  gtc <- hcocena::GeneToCluster()
  genes <- dplyr::filter(gtc, color == cluster) %>% 
    dplyr::pull(., "gene")
  return(genes)
}


#' Internal Function Used In zoom()
#' @noRd

filter_counts <- function(genes){
  out <- list()
  for(x in base::seq_along(hcobject[["layers"]])){
    out[[base::paste0("set", x, "_counts")]] <- dplyr::filter(hcobject[["data"]][[base::paste0("set", x, "_counts")]], rownames(hcobject[["data"]][[base::paste0("set", x, "_counts")]]) %in% genes)
  }
  return(out)
}



#' Internal Function Used In zoom()
#' @noRd

filter_correlations <- function(genes){
  out <- list()
  for(x in base::seq_along(hcobject[["layers"]])){
    tmp <- hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["corr_calc_out"]][["correlation_df_filt"]]
    out[[base::paste0("corrdf_", x)]] <- dplyr::filter(tmp, V1 %in% genes & V2 %in% genes)
  }
  return(out)
}



#' Internal Function Used In zoom()
#' @noRd

get_cutoff_stats <- function(corr_dfs, min_nodes){
  out <- list()
  for(x in base::seq_along(hcobject[["layers"]])){
    cutoff_stats <- base::do.call("rbind", base::lapply(X = hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["corr_calc_out"]][["range_cutoff"]],
                                                        FUN = cutoff_prep,
                                                        corrdf_r = corr_dfs[[base::paste0("corrdf_", x)]],
                                                        print.all.plots = FALSE,
                                                        x = x,
                                                        min_nodes = min_nodes))
    # reshape cutoff stats:
    out[[paste0("cutoff_calc_out_", x)]] <- reshape_cutoff_stats(cutoff_stats = cutoff_stats)
  }
  return(out)
}



#' Internal Function Used In zoom()
#' @noRd

plot_cutoffs_zoom <- function(interactive = T, 
                             hline = list("R.squared" = NULL, 
                                          "no_edges" = NULL, 
                                          "no_nodes" = NULL, 
                                          "no_networks" = NULL),
                             cutoff_stats){
  
  for(x in 1:base::length(hcobject[["layers"]])){
    cutoff_df <- cutoff_stats[[paste0("cutoff_calc_out_", x)]][["cutoff_stats_concise"]]
    if(interactive == T){
      p <- plot_cutoffs_internal_interactive(cutoff_stats = cutoff_df,
                                             x = x)
    }else{
      p <- plot_cutoffs_internal_static(cutoff_stats = cutoff_df, 
                                        hline = hline, 
                                        x = x)
    }
    
    print(p)
    
  }
}	


#' Zoom Into CLuster Part I
#' 
#' Reruns all CoCena steps until the cutoff selection for a selected cluster.
#' @export

zoom <- function(){
  out <- list()

  cluster <- hcobject[["satellite_outputs"]][["zoom"]][["settings"]][["cluster"]]
  message("... Zooming  into cluster ", cluster, " ...")
  out[["genes"]] <- get_cluster_genes(cluster = cluster)
  out[["filt_counts"]] <- filter_counts(genes = out[["genes"]])
  out[["filt_corrs"]] <- filter_correlations(genes = out[["genes"]])
  out[["cutoff_stats"]] <- get_cutoff_stats(corr_dfs = out[["filt_corrs"]], min_nodes = hcobject[["satellite_outputs"]][["zoom"]][["settings"]][["min_nodes"]])
  plot_cutoffs_zoom(interactive = T, cutoff_stats = out[["cutoff_stats"]])
  hcobject[["satellite_outputs"]][["zoom"]][[cluster]] <<- out
}



#' Internal Function Used In zoom2()
#' @noRd

zoom_network_construction <- function(zoom1res, cutoff_vec = rep(0, base::length(hcobject[["layers"]])), min_nodes){
  out <- list()
  for(x in base::seq_along(hcobject[["layers"]])){
    
    # extract data exceeding set cutoff:
    filt_cutoff_data <- zoom1res[["filt_corrs"]][[paste0("corrdf_", x)]] %>% 
      dplyr::filter(., rval >= cutoff_vec[x])
    
    # build network based on filtered data:
    filt_cutoff_graph <- igraph::graph_from_data_frame(filt_cutoff_data, directed = FALSE)
    
    # remove too small components:
    graph_components <- igraph::components(filt_cutoff_graph)
    
    # remove nodes from too small components from the network:
    gene_to_comp <- base::data.frame(gene = base::names(graph_components[["membership"]]), component = graph_components[["membership"]])
    comps_to_keep <- base::which(graph_components[["csize"]] >= min_nodes)
    nodes_to_remove <- dplyr::filter(gene_to_comp, !component %in% comps_to_keep) %>% dplyr::pull(., "gene")
    filt_cutoff_graph <- igraph::delete.vertices(filt_cutoff_graph, nodes_to_remove)
    
    # remove edges from removed nodes from the cutoff data:
    filt_cutoff_data <- dplyr::filter(filt_cutoff_data, V1 %in% igraph::V(filt_cutoff_graph)$name & V2 %in% igraph::V(filt_cutoff_graph)$name)
    
    out[[paste0("filt_cutoff_graph_", x)]] <- filt_cutoff_graph
    out[[paste0("filt_cutoff_data_", x)]] <- filt_cutoff_data
    
  }
  return(out)
}



#' Internal Function Used In zoom2()
#' @noRd

integrate_zoom_network <- function(zn, mode = 'u', multi_edges = "min", with = NULL){
  
  out <- list()
  
  if(mode == "u"){
    message("Intergrating network based on union.")
    combined_edgelist <- NULL
    
    for (x in base::seq_along(hcobject[["layers"]])){
      combined_edgelist <- base::rbind(combined_edgelist, zn[[paste0("filt_cutoff_data_", x)]])
    }
    combined_edgelist$weight <- combined_edgelist$rval
    combined_edgelist$rval <- NULL
    combined_edgelist$pval <- NULL
    
  }else if(mode == "i"){
    message("Intergrating network based on intersection.")
    if(is.null(with)){
      stop("The 'with' parameter must be specified.")
    }
    # change to numeric representation of the dataset in case it was given as the name of a dataset:
    if(with %in% hcobject[["layers_names"]]){
      with <- match(with, hcobject[["layers_names"]])
    }
    
    # change to numeric representation of the dataset in case it was given as a string:
    if(startsWith(with, "set")){
      with <- base::strsplit(base::as.character(with), split = "set")[[1]][2] %>% base::as.numeric()
    }
    
    # get edgelist of the reference network:
    combined_edgelist <- zn[[paste0("filt_cutoff_data_", with)]]
    
    combined_edgelist$merged <- base::paste0(combined_edgelist$V1 %>% base::as.character(), 
                                             combined_edgelist$V2 %>% base::as.character())
    combined_edgelist$revmerged <- base::paste0(combined_edgelist$V2 %>% base::as.character(), 
                                                combined_edgelist$V1 %>% base::as.character())
    # iterate over datasets:
    for (x in base::seq_along(hcobject[["layers"]])){
      
      if(!x == with){
        # get edgelist of current dataset:
        tmp <- zn[[paste0("filt_cutoff_data_", x)]]
        
        tmp$merged <- base::paste0(tmp$V1 %>% base::as.character(), 
                                   tmp$V2 %>% base::as.character())
        
        # check which edges overlap with reference network:
        tmp <- dplyr::filter(tmp, merged %in% combined_edgelist$merged | merged %in% combined_edgelist$revmerged)
        combined_edgelist <- base::rbind(combined_edgelist, tmp)
      }
    }
    combined_edgelist$weight <- combined_edgelist$rval
    combined_edgelist$rval <- NULL
    combined_edgelist$pval <- NULL
    combined_edgelist$merged <- NULL
    combined_edgelist$revmerged <- NULL
    
  }else{
    stop("No valid choice of 'mode'. Must be either 'u' for integration by union or 'i' for integration by intersection.")
  }
  merged_net <- igraph::graph_from_data_frame(combined_edgelist, directed=FALSE)
  merged_net <- igraph::simplify(merged_net, edge.attr.comb=list(weight=multi_edges, "ignore"))
  new_edgelist <- base::cbind( igraph::get.edgelist(merged_net) , base::round(igraph::E(merged_net)$weight, 7)) %>% base::as.data.frame()
  base::colnames(new_edgelist) <- base::colnames(combined_edgelist)
  out[["combined_edgelist"]] <- new_edgelist
  out[["merged_net"]] <- merged_net
  return(out)
}



#' Internal Function Used In zoom2()
#' @noRd

cluster_calculation_zoom <- function(zint,
                                cluster_algo,
                                no_of_iterations = 1,
                                max_cluster_count_per_gene = 1,
                                min_nodes_number_for_cluster = 5) {
  
  
  if(cluster_algo == "cluster_leiden"){
    cluster_information <- leiden_clustering_zoom(g = zint[["merged_net"]], num_it = no_of_iterations, min_nodes_number_for_cluster = min_nodes_number_for_cluster)
    cluster_information <- dplyr::filter(cluster_information, !color == "white")
    return(cluster_information)
  }else{
    
    # count the number of graph components present in the network:
    comps <- igraph::count_components(zint[["merged_net"]])
    
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
      df_modularity_score <- base::do.call("rbind", base::lapply(algos_to_use,
                                                                 cluster_calculation_internal_zoom,
                                                                 graph_obj = zint[["merged_net"]],
                                                                 case = "test",
                                                                 min_nodes_number_for_cluster = min_nodes_number_for_cluster))
      
      
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
    igraph_list[[1]] <- zint[["merged_net"]]
    
    # apply best clustering algorithm:
    gene_which_cluster <- base::do.call("cbind", base::lapply(1:no_of_iterations,
                                                              cluster_calculation_internal_zoom,
                                                              algo = cluster_algo_used,
                                                              case = "best",
                                                              graph_obj = zint[["merged_net"]]))
    
    
    
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
    
    
    
    cluster_Data <- base::data.frame(genes = igraph::vertex_attr(zint[["merged_net"]], "name"),
                                     clusters = base::paste0("Cluster ", gene_cluster_ident),
                                     stringsAsFactors = FALSE)
    
    #summarize the data
    #produces a table where col are cluster name, number of components,
    #names of genes in cluster
    
    dfk <- cluster_Data %>%
      dplyr::count(clusters,genes) %>%
      dplyr::group_by(clusters) %>%
      dplyr::summarise(gene_no = base::sum(n), gene_n = base::paste0(genes, collapse = ",")) %>%
      dplyr::mutate(cluster_included = base::ifelse(gene_no >= min_nodes_number_for_cluster, "yes", "no"), color = "white")
    
    
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
                       min_nodes_number_for_cluster, " genes per cluster. These genes will also be assigned to Cluster 0 (white) and left out of the network and further analyses."))
    
    
    
    
    # for each cluster produces a row of means per condition (info data) for all genes within the cluster
    
    gfc_dat <- hcobject[["integrated_output"]][["GFC_all_layers"]]
    
    dfk_allinfo <- base::do.call("rbind", base::lapply(1:base::nrow(dfk), gfc_mean_clustergene, 
                                                       cluster_df = dfk,
                                                       gfc_dat = gfc_dat))
    
    dfk_allinfo$vertexsize <- base::ifelse(dfk_allinfo$cluster_included == "yes", 3, 1)

    dfk_allinfo <- dplyr::filter(dfk_allinfo, !color == "white")

    
    return(dfk_allinfo)
    
  }
  
}



#' Internal Function Used In zoom2()
#' @noRd

leiden_clustering_zoom <- function(g, num_it, min_nodes_number_for_cluster){
  
  color.cluster <- get_cluster_colours()
  
  set.seed(168575)
  
  # run Leiden algorithm ion network:
  partition <- leidenAlg::leiden.community(graph = g, n.iterations = num_it)
  
  # extract found clusters and the genes belonging to them:
  clusters_df <- base::data.frame(cluster = base::as.numeric(partition$membership), gene = partition$names)
  
  # the algo start enumeration at 0, increase to 1 for easier comaptibility with R fucntions:
  clusters_df$cluster <- clusters_df$cluster + 1
  
  # get gene counts per cluster:
  cluster_frequencies <- base::table(clusters_df$cluster) %>% base::as.data.frame() 
  
  # detect clusters large enough to be kept:
  clusters_to_keep <- dplyr::filter(cluster_frequencies, Freq >= min_nodes_number_for_cluster) %>% 
    dplyr::pull(., "Var1")
  
  # define white clusters (those that are too small to be kept):
  clusters_df_white <- dplyr::filter(clusters_df, !cluster %in% clusters_to_keep)
  
  # remove white clusters:
  clusters_df <- dplyr::filter(clusters_df, cluster %in% clusters_to_keep)
  
  # inform how many clusters and accordingly how many genes were lost due to insufficient cluster size:
  print(base::paste0(base::length(base::unique(clusters_df_white$cluster)), 
                     " cluster/s was/were smaller than the set minimum cluster size and therefore discarded.",
                     " This removes ", base::nrow(clusters_df_white), " genes."))
  
  
  out <- base::lapply(base::unique(clusters_df$cluster), function(x){
    
    # extract genes present in current cluster:
    genes <- dplyr::filter(clusters_df, cluster == x) %>% 
      dplyr::pull(., "gene")
    
    # cluster name:
    col_clusters <- base::paste0("cluster ", x)
    # number of genes in cluster:
    col_gene_no <- dplyr::filter(clusters_df, cluster == x) %>% 
      base::nrow()
    # comma separated list of all genes in the cluster:
    col_gene_n <- genes %>% base::paste0(., collapse = ",")
    # is the cluster included in the network (i.e., is it large enough):
    col_cluster_included <- "yes"
    # cluster color:
    col_color <- color.cluster[x]
    # order of conditions for following GFC values:
    col_conditions <- base::colnames(dplyr::select(hcobject[["integrated_output"]][["GFC_all_layers"]], -Gene)) %>% 
      base::paste0(., collapse = "#")
    # mean GFCs of the cluster genes per sample group:
    col_grp_means <- dplyr::filter(hcobject[["integrated_output"]][["GFC_all_layers"]], Gene %in% genes) %>% 
      dplyr::select(-Gene) %>% 
      base::apply(., 2, base::mean) %>% 
      base::round(., digits = 3) %>% 
      base::paste0(., collapse = ",")
    # collect all information:
    out <- base::data.frame(clusters = col_clusters,
                            gene_no = col_gene_no,
                            gene_n = col_gene_n,
                            cluster_included = col_cluster_included,
                            color = col_color,
                            conditions = col_conditions,
                            grp_means = col_grp_means,
                            vertexsize = 3)
    return(out)
  }) %>% rlist::list.rbind()
  return(out)
}




#' Internal Function Used In zoom2()
#' @noRd

cluster_calculation_internal_zoom <- function(graph_obj, 
                                              algo, 
                                              case, 
                                              it = 1,
                                              min_nodes_number_for_cluster) {
  
  library(igraph)
  
  if(algo == "cluster_leiden"){
    cfg <- leiden_clustering_zoom(g = graph_obj, num_it = no_of_iterations, min_nodes_number_for_cluster = min_nodes_number_for_cluster)
  }
  cfg <- base::get(algo)(graph_obj)
  
  mod_score <- igraph::modularity(graph_obj, base::as.numeric(cfg$membership)) 
  
  mod_df <- base::data.frame(modularity_score = mod_score, 
                             cluster_algorithm = algo, 
                             stringsAsFactors = F)
  
  # making switch so that in the end when only the best algorithm is to be used then the same function can be used
  output <- base::switch(case, best = cfg$membership, test = mod_df, final = cfg)
  
  print(base::paste0(algo," algorithm tested"))
  return(output)
}




#' Internal Function Used In zoom2()
#' @noRd

plot_cluster_heatmap_zoom <- function(cluster_information,
                                 col_order = NULL, 
                                 row_order = NULL, 
                                 cluster_columns = T,
                                 cluster_rows = T, 
                                 k = 0,  
                                 cat_as_bp = NULL, 
                                 file_name = "zoom_module_heatmap.pdf"){
  
 
  
  # get categorical and numerical sample group annotations if they exist:
  if("column_annos_categorical" %in% base::names(hcobject[["satellite_outputs"]])){
    column_anno_categorical <- hcobject[["satellite_outputs"]][["column_annos_categorical"]]
  }else{
    column_anno_categorical <- NULL
  }
  
  if("column_annos_numerical" %in% base::names(hcobject[["satellite_outputs"]])){
    column_anno_numerical <- hcobject[["satellite_outputs"]][["column_annos_numerical"]]
  }else{
    column_anno_numerical <- NULL
  }
  
  if(base::is.null(cat_as_bp)){
    if(!base::is.null(column_anno_categorical)){
      cat_as_bp <- base::rep(F, base::length(column_anno_categorical))
    }
  }
  
  base::gc()
  
  # filter for included clusters (non-white)
  c_df <- dplyr::filter(cluster_information, cluster_included == "yes")
  mat_heatmap <- NULL
  
  
  if(!base::is.null(row_order)){
    for (c in row_order){
      #get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)
      
      # GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(hcobject[["integrated_output"]][["GFC_all_layers"]], Gene %in% genes)
      c_GFC_means <- base::apply(c_GFCs[, base::c(1:(base::ncol(c_GFCs)-1))], 2, base::mean)
      
      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)
      
    }
    base::rownames(mat_heatmap) <- row_order
    
  }else{
    for (c in base::unique(c_df$color)){
      #get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)
      
      
      # GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(hcobject[["integrated_output"]][["GFC_all_layers"]], Gene %in% genes)
      
      if(base::is.vector(c_GFCs)){
        c_GFC_means <- cGFCs
      }else{
        c_GFC_means <- base::apply(c_GFCs[, base::c(1:(base::ncol(c_GFCs)-1))] %>% 
                                     base::as.data.frame(), 2, base::mean)
      }
      
      
      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)
      
    }
    base::rownames(mat_heatmap) <- c_df$color
  }
  
  
  base::colnames(mat_heatmap) <- base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])[1:(base::ncol(hcobject[["integrated_output"]][["GFC_all_layers"]])-1)]
  
  if(!base::is.null(col_order)){
    mat_heatmap <- mat_heatmap %>% base::as.data.frame()
    mat_heatmap <- mat_heatmap[, col_order] %>% base::as.matrix()
  }
  
  enrich_mat1 <- list()
  enrich_count1 <- list()
  enrich_mat2 <- list()
  enrich_count2 <- list()
  
  
  
  if(!base::length(enrich_mat1) == 0){
    enrich_mat1 <- base::matrix(base::unlist(enrich_mat1), nrow = base::length(enrich_mat1), byrow = T)
    enrich_count1 <- base::unlist(enrich_count1)
    
  }
  if(base::all(enrich_mat1 == 0) == T){
    enrich_mat1 <- list()
  }
  
  if(!base::length(enrich_mat2) == 0){
    
    enrich_mat2 <- base::matrix(base::unlist(enrich_mat2), nrow = base::length(enrich_mat2), byrow = T)
    enrich_count2 <- base::unlist(enrich_count2)
  }
  if(base::all(enrich_mat2 == 0) == T){
    enrich_mat2 <- list()
  }
  
  if(!base::is.null(row_order)){
    cluster_colors <- base::factor(row_order)
    base::names(cluster_colors) <- row_order
    c_df <- c_df[base::match(row_order, c_df$color),]
  }else{
    cluster_colors <- base::factor(c_df$color)
    base::names(cluster_colors) <- c_df$color
    row_order <- base::unique(c_df$color)
  }
  
  
  if(base::length(enrich_mat1) == 0 & base::length(enrich_mat2) == 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = ComplexHeatmap::anno_simple(row_order, col = cluster_colors,
                                                                                  simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")),
                                            genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(2.5, "cm")),
                                            gene_nums = ComplexHeatmap::anno_text(c_df$gene_no, width = grid::unit(1.5, "cm"), gp = grid::gpar(fontsize = 10)),
                                            
                                            
                                            which = "row",
                                            width = grid::unit(4.5, "cm"),
                                            annotation_name_side = "top",
                                            gap = grid::unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = grid::gpar(fontsize = 8))
    
    lgd_list <- list(
      
    )
  }
  else if(base::length(enrich_mat1) > 0 & base::length(enrich_mat2) == 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = ComplexHeatmap::anno_simple(row_order, col = cluster_colors,
                                                                                  simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")),
                                            genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(2.5, "cm")),
                                            
                                            enriched_count = ComplexHeatmap::anno_text(base::paste0(enrich_count1, "/", c_df$gene_no), width = grid::unit(1.5, "cm")),
                                            enriched = ComplexHeatmap::anno_barplot(enrich_mat1,
                                                                                    width = grid::unit(3, "cm"),
                                                                                    gp = grid::gpar(fill = RColorBrewer::brewer.pal(n = 12, name = "Paired"),
                                                                                                    col = RColorBrewer::brewer.pal(n = 12, name = "Paired"))),
                                            which = "row",
                                            width = grid::unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = grid::unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = grid::gpar(fontsize = 8))
    
    lgd_list <- list(
      
      ComplexHeatmap::Legend(labels = base::unique(user_enrichment_1$cell_type), title = "enriched",
                             legend_gp = grid::gpar(col = RColorBrewer::brewer.pal(n = 12, name = "Paired")),
                             type = "points", pch = 15)
    )
  }
  else if(base::length(enrich_mat1) == 0 & base::length(enrich_mat2) > 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = ComplexHeatmap::anno_simple(row_order, col = cluster_colors,
                                                                                  simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")),
                                            genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(1.5, "cm")),
                                            
                                            enriched_count = ComplexHeatmap::anno_text(base::paste0(enrich_count2, "/", c_df$gene_no), width = grid::unit(1.5, "cm"),
                                                                                       gp = grid::gpar(fontsize = 8)),
                                            enriched = ComplexHeatmap::anno_barplot(enrich_mat2,
                                                                                    width = grid::unit(5, "cm"),
                                                                                    gp = grid::gpar(fill = ggsci::pal_d3(palette = "category20")(base::ncol(enrich_mat2)),
                                                                                                    col = ggsci::pal_d3(palette = "category20")(base::ncol(enrich_mat2)))),
                                            
                                            which = "row",
                                            width = grid::unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = grid::unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = grid::gpar(fontsize = 8))
    
    lgd_list <- list(
      
      
      ComplexHeatmap::Legend(labels = base::unique(user_enrichment_2$cell_type), title = "enriched",
                             legend_gp = grid::gpar(col = ggsci::pal_d3(palette = "category20")(20)),
                             type = "points", pch = 15)
    )
  }else{
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = ComplexHeatmap::anno_simple(row_order, col = cluster_colors,
                                                                                  simple_anno_size = grid::unit(0.25, "cm"), gp = grid::gpar(col = "black")),
                                            genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(0.75, "cm")),
                                            
                                            enriched_count_1 = ComplexHeatmap::anno_text(base::paste0(enrich_count1, "/", c_df$gene_no),
                                                                                         width = grid::unit(0.75, "cm"),
                                                                                         gp = grid::gpar(fontsize = 8)),
                                            enriched_1 = ComplexHeatmap::anno_barplot(enrich_mat1,
                                                                                      width = grid::unit(2, "cm"),
                                                                                      gp = grid::gpar(fill = RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:base::ncol(enrich_mat1)],
                                                                                                      col = RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:base::ncol(enrich_mat1)]),
                                                                                      baseline = 0),
                                            enriched_count_2 = ComplexHeatmap::anno_text(base::paste0(enrich_count2, "/", c_df$gene_no),
                                                                                         width = grid::unit(0.75, "cm"),
                                                                                         gp = grid::gpar(fontsize = 8)),
                                            enriched_2 = ComplexHeatmap::anno_barplot(enrich_mat2,
                                                                                      width = grid::unit(2, "cm"),
                                                                                      gp = grid::gpar(fill = ggsci::pal_d3(palette = "category20")(20)[1:base::ncol(enrich_mat2)],
                                                                                                      col = ggsci::pal_d3(palette = "category20")(20)[1:base::ncol(enrich_mat2)]),
                                                                                      baseline = 0),
                                            which = "row",
                                            width = grid::unit(12, "cm"),
                                            annotation_name_side = "top",
                                            gap = grid::unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = grid::gpar(fontsize = 8))
    
    lgd_list <- list(
      
      
      ComplexHeatmap::Legend(labels = base::unique(user_enrichment_1$cell_type), title = "enriched_1",
                             legend_gp = grid::gpar(col = ggsci::pal_d3(palette = "category20")(20)[1:base::ncol(enrich_mat1)]),
                             type = "points", pch = 15),
      ComplexHeatmap::Legend(labels = base::unique(user_enrichment_2$cell_type), title = "enriched_2",
                             legend_gp = grid::gpar(col = ggsci::pal_d3(palette = "category20")(20)[1:base::ncol(enrich_mat2)]),
                             type = "points", pch = 15)
    )
  }
  
  
  anno_list <- NULL
  
  
  if(!base::length(column_anno_categorical) == 0){
    for(a in 1:base::length(column_anno_categorical)){
      base::set.seed(a)
      tmp_colour <- grDevices::colorRampPalette(ggsci::pal_d3("category20")(20))(base::ncol(column_anno_categorical[[a]]))
      if(cat_as_bp[a] == T){
        column_anno_categorical[[a]][base::is.na(column_anno_categorical[[a]])] <- 0
        if(base::is.null(anno_list)){
          anno_list <- ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_barplot(column_anno_categorical[[a]]%>%base::as.matrix(),
                                                                                                 width = grid::unit(2, "cm"),
                                                                                                 gp = grid::gpar(fill = tmp_colour,
                                                                                                                 col = tmp_colour)),
                                                         which = "column",
                                                         height = grid::unit(1, "cm"),
                                                         annotation_name_side = "right",
                                                         gap = grid::unit(2, "mm"),
                                                         annotation_name_rot = 0,
                                                         annotation_name_gp = grid::gpar(fontsize = 8),
                                                         annotation_label = base::names(column_anno_categorical)[a])
        }else{
          anno_list <-  ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_barplot(column_anno_categorical[[a]]%>%base::as.matrix(),
                                                                                                                                         width = grid::unit(2, "cm"),
                                                                                                                                         gp = grid::gpar(fill = tmp_colour,
                                                                                                                                                         col = tmp_colour)),
                                                                                                 which = "column",
                                                                                                 height = grid::unit(1, "cm"),
                                                                                                 annotation_name_side = "right",
                                                                                                 gap = grid::unit(2, "mm"),
                                                                                                 annotation_name_rot = 0,
                                                                                                 annotation_name_gp = grid::gpar(fontsize = 8),
                                                                                                 annotation_label = base::names(column_anno_categorical)[a]), direction = c("vertical"))
        }
        
        
        
      }else{
        if(base::is.null(anno_list)){
          anno_list <-ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_lines(column_anno_categorical[[a]] %>% base::as.matrix(),width = grid::unit(2, "cm"),
                                                                                              gp = grid::gpar(col = tmp_colour),
                                                                                              add_points = TRUE,
                                                                                              pt_gp = grid::gpar(col = tmp_colour), pch = 16),
                                                        which = "column",
                                                        height = grid::unit(1, "cm"),
                                                        annotation_name_side = "right",
                                                        gap = grid::unit(2, "mm"),
                                                        annotation_name_rot = 0,
                                                        annotation_name_gp = grid::gpar(fontsize = 8),
                                                        annotation_label = base::names(column_anno_categorical)[a])
          
        }else{
          anno_list <-  ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_lines(column_anno_categorical[[a]] %>% base::as.matrix(),width = grid::unit(2, "cm"),
                                                                                                                                       gp = grid::gpar(col = tmp_colour),
                                                                                                                                       add_points = TRUE,
                                                                                                                                       pt_gp = grid::gpar(col = tmp_colour), pch = 16),
                                                                                                 which = "column",
                                                                                                 height = grid::unit(1, "cm"),
                                                                                                 annotation_name_side = "right",
                                                                                                 gap = grid::unit(2, "mm"),
                                                                                                 annotation_name_rot = 0,
                                                                                                 annotation_name_gp = grid::gpar(fontsize = 8),
                                                                                                 annotation_label = base::names(column_anno_categorical)[a]), direction = c("vertical"))
        }
        
      }
      
      
      
      lgd_list <- rlist::list.append(lgd_list, ComplexHeatmap::Legend(labels = base::colnames(column_anno_categorical[[a]]%>%base::as.matrix()), title = base::names(column_anno_categorical)[a],
                                                                      legend_gp = grid::gpar(col = tmp_colour),
                                                                      type = "points", pch = 15))
    }
    
  }
  
  
  
  if(!base::length(column_anno_numerical) == 0){
    for(a in 1:base::length(column_anno_numerical)){
      tmp_col_anno_2 <- column_anno_numerical[[a]]
      tmp_col_anno_2 <- tmp_col_anno_2[base::colnames(mat_heatmap)]
      if(base::is.null(anno_list)){
        anno_list <- ComplexHeatmap::HeatmapAnnotation(cont_anno = ComplexHeatmap::anno_boxplot(tmp_col_anno_2, height = grid::unit(1, "cm")),
                                                       which = "column",
                                                       annotation_name_side = "right",
                                                       gap = grid::unit(2, "mm"),
                                                       annotation_name_rot = 0,
                                                       annotation_name_gp = grid::gpar(fontsize = 8),
                                                       annotation_label = base::names(column_anno_numerical)[a], show_legend = F)
        
      }else{
        anno_list <-  ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::HeatmapAnnotation(cont_anno = ComplexHeatmap::anno_boxplot(tmp_col_anno_2, height = grid::unit(1, "cm")),
                                                                                               which = "column",
                                                                                               annotation_name_side = "right",
                                                                                               gap = grid::unit(2, "mm"),
                                                                                               annotation_name_rot = 0,
                                                                                               annotation_name_gp = grid::gpar(fontsize = 8),
                                                                                               annotation_label = base::names(column_anno_numerical)[a], show_legend = F), direction = c("vertical"))
      }
      
      
    }
  }
  
  
  all_conditions <- NULL
  
  # if(!is.null(column_anno_categorical) | !is.null(column_anno_numerical)){
  for(setnum in 1:base::length(hcobject[["layers"]])){
    all_conditions <- base::c(all_conditions, base::as.character(dplyr::pull(hcobject[["data"]][[base::paste0("set", setnum, "_anno")]], hcobject[["global_settings"]][["voi"]])))
  }
  all_conditions <- base::table(all_conditions) %>%
    base::as.data.frame() %>%
    dplyr::filter(., all_conditions %in% base::colnames(mat_heatmap))
  all_conditions <- all_conditions[base::match(base::colnames(mat_heatmap), base::as.character(all_conditions$all_conditions)),]
  all_conditions <- base::paste0(all_conditions$all_conditions, "  [", all_conditions$Freq, "]")
  if(base::is.null(anno_list)){
    anno_list <- ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions))
    
  }else{
    anno_list <-  ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions)), direction = c("vertical"))
  }
  
  
  # }
  
  
  
  Cairo::Cairo(file = paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", file_name),
               width = 50,
               height = 30,
               pointsize=11,
               dpi=300,
               type = "pdf",
               units = "in")
  
  hm <- ComplexHeatmap::Heatmap(mat_heatmap,
                                right_annotation = ha,
                                col = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(base::length(base::seq(-2, 2, by = .1))),
                                clustering_distance_rows = "euclidean",
                                clustering_distance_columns = "euclidean",
                                clustering_method_rows = "complete",
                                clustering_method_columns = "complete",
                                cluster_columns = cluster_columns,
                                cluster_rows = cluster_rows,
                                column_names_rot = 90,
                                column_names_centered = F,
                                row_names_gp = grid::gpar(fontsize = 8),
                                column_names_gp = grid::gpar(fontsize = 8),
                                rect_gp = grid::gpar(col = "black"),
                                heatmap_legend_param = list(title = "", legend_height = grid::unit(3, "cm")), column_km = k)
  if(base::is.null(anno_list)){
    anno_list <- hm 
  }else{
    anno_list <- ComplexHeatmap::add_heatmap(hm, anno_list, direction = c("vertical"))
  }
  
  hm_w_lgd <- ComplexHeatmap::draw(anno_list, annotation_legend_list = lgd_list, merge_legends = T,
                                   padding = grid::unit(c(2, 2, 2, 30), "mm"))
  
  grDevices::dev.off()
  
  print(hm_w_lgd)

  return(hm_w_lgd)
}




#' Internal Function Used In zoom2()
#' @noRd

plot_zoom_network <- function(network,
                              cluster_info,
                              layout = NULL, 
                              gene_labels = NULL, 
                              save = F, 
                              label_offset = 50){
  

  gene_to_cluster <- GeneToCluster(cluster_information = cluster_info)
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
    message("... Calculating Fruchtermann Rheingold Layout ...")
    l <- igraph::layout.fruchterman.reingold(network)
    base::rownames(l) <- igraph::V(network)$name
    
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
      Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/zoom_network_col_by_module.pdf"), 
                      width = 15, height = 15)
      
      igraph::plot.igraph(network2, vertex.size = 3, vertex.label = new_labels, vertex.label.cex = 0.75,
                          layout = l2, main = "zoom network coloured by module",
                          edge.color = new_edge_color, vertex.label.color = new_label_color,
                          vertex.frame.color = new_frame_color)
      
      grDevices::dev.off()
    }
    
    
    igraph::plot.igraph(network2, vertex.size = 3, vertex.label = new_labels, vertex.label.cex = 0.75,
                        layout = l2, main = "zoomed network coloured by module",
                        edge.color = new_edge_color, vertex.label.color = new_label_color,
                        vertex.frame.color = new_frame_color)
    
  }else{
    
    if(save == T){
      Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/zoom_network_col_by_module.pdf"), 
                      width = 15, height = 15)
      
      igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3, 
                          layout = l, main = "zoomed network coloured by cluster")
      
      grDevices::dev.off()
    }
    
    
    igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3, 
                        layout = l, main = "zoomed network coloured by cluster")
    
  }

  
}




#' Cluster Zoom Settings
#'  
#' Sets the parameters to zoom into a specific cluster.
#' @param cluster The colour of the cluster as a string.
#' @param min_nodes An integer. The minimum number of nodes to comprise a cluster, default is 15.
#' @param mode How to integrate the networks of the different layers (can be ignored, if there is only one dataset). Default is 'u' (integration by union), alternatively you can set 'i' (integration by intersection).
#' @param multi_edges One of “min” (default), “mean” or “max” resulting in the simplification of the multigraph by using the minimum, the mean or the maximum edge weight among the multiple edges, respectively.
#'  Multiple edges occur when the edge was present in more than one datset.
#' @param no_of_iterations How often to repeat the clsutering algorithm. Default is 10.
#' @param max_cluster_count_per_gene The maximum number of clsuters a gene can be associated with over the different iteration before it is discared into cluster 'white'. Default is 1.
#' @param min_nodes_number_for_cluster The minimum number of nodes to comprise a cluster in the network. Default is 5.
#' @export

zoom_settings <- function(cluster,
                          min_nodes = 15,
                          mode = "u",
                          multi_edges = "min",
                          no_of_iterations = 10,
                          max_cluster_count_per_gene = 1,
                          min_nodes_number_for_cluster = 5){

  hcobject[["satellite_outputs"]][["zoom"]][["settings"]] <<- list(cluster = cluster,
                          min_nodes = min_nodes,
                          mode = mode,
                          multi_edges = multi_edges,
                          no_of_iterations = no_of_iterations,
                          max_cluster_count_per_gene = max_cluster_count_per_gene,
                          min_nodes_number_for_cluster = min_nodes_number_for_cluster)

}




#' Zoom Into CLuster Part II
#' 
#' Filters Cluster genes for set cutoff, integrates layers (if necessary), re-clusters and re-plots cluster-network
#' @param cutoff_vec A vector of cutoff values. The vector's length must be equal to the number of datasets. The order in which the cutoffs are set must correspond to the order in which the layers were declared in define_layers().
#' @param cluster_algo The clustering algorithm to be used. The choise is between "cluster_leiden" (default), "cluster_louvain", "cluster_label_prop", "cluster_fast_greedy",
#'  "cluster_infomap", "cluster_walktrap" and "auto" (in which case all are tested and the one with the highest modularity is chosen).
#' @export

zoom2 <- function(cutoff_vec, cluster_algo = "cluster_leiden"){

  cluster <- hcobject[["satellite_outputs"]][["zoom"]][["settings"]][["cluster"]]
  
  zn <- zoom_network_construction(zoom1res = hcobject[["satellite_outputs"]][["zoom"]][[cluster]], cutoff_vec = cutoff_vec, min_nodes = hcobject[["satellite_outputs"]][["zoom"]][["settings"]][["min_nodes"]])
  
  
  zint <- integrate_zoom_network(zn = zn, mode = hcobject[["satellite_outputs"]][["zoom"]][["settings"]][["mode"]], multi_edges = hcobject[["satellite_outputs"]][["zoom"]][["settings"]][["multi_edges"]])
  
  
  zclustinfo <- cluster_calculation_zoom(zint = zint, 
                                         cluster_algo = cluster_algo,
                                         no_of_iterations = hcobject[["satellite_outputs"]][["zoom"]][["settings"]][["no_of_iterations"]], 
                                         max_cluster_count_per_gene = hcobject[["satellite_outputs"]][["zoom"]][["settings"]][["max_cluster_count_per_gene"]], 
                                         min_nodes_number_for_cluster = hcobject[["satellite_outputs"]][["zoom"]][["settings"]][["min_nodes_number_for_cluster"]])
  
  
  zhm <- plot_cluster_heatmap_zoom(cluster_information = zclustinfo)
  
  
  plot_zoom_network(network = zint$merged_net, cluster_info = zclustinfo)
  
  out <- list()
  hcobject[["satellite_outputs"]][["zoom"]][[cluster]][["clust_info"]] <<- zclustinfo
  hcobject[["satellite_outputs"]][["zoom"]][[cluster]][["network"]] <<- zint$merged_net
  hcobject[["satellite_outputs"]][["zoom"]][[cluster]][["heatmap"]] <<- zhm
}