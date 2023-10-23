#' Compare PCAs For Different Clustering Algorithms
#' 
#' A series of principle component analyses (PCAs) are performed and plotted. The first one to be shown is based on the top most variant genes (or all depending on the topvar setting) for each data set. 
#' 	The following PCAs demonstrate for each clustering algorithm the spatial arrangement of the samples in the space of the first two principle components based on their expression behaviour among the clusters. 
#' 	Samples that show similar behaviour among the defined clusters are located closer together than samples that show very different expressions across clusters. 
#' 	A suitable clustering algorithm will identify clusters such that samples that have similar underlying data (first PCA) will also be closely situated in the PCA based on similar cluster expressions.
#' @param algo If NULL (default), PCAs for all available algorithms are plotted. 
#' 	Otherwise, one of "cluster_louvain", "cluster_fast_greedy", "cluster_infomap", "cluster_walktrap", "cluster_label_prop", "cluster_leiden", then only the PCA for that algorithm is shown.
#' @param gtc A user defined clsutering can be provided as long as it has the same format as the output of GeneToCluster(): A data frame with two columns, the first containing gene names as strings, the second containing cluster colours as strings.
#' @param cols A user-defined color vector.
#' @export

PCA_algo_compare <- function(gtc = NULL, algo = NULL, cols = NULL){

  # recycle clustering of all algos from alluvial function, if it exists, otherwise conduct clustering now:
  if(!"alluvials" %in% names(hcobject[["integrated_output"]])){
  	all_clusterings <- run_all_cluster_algos()
  }else{
  	all_clusterings <- hcobject[["integrated_output"]][["alluvials"]]
  }

  PCAs <- list()
  PCAs[["topvar"]] <- plot_PCA_topvar(PCA_save_folder = "pca_algo_compare", cols = cols)

  if(base::is.null(gtc)){
    if(base::is.null(algo)){
      for(a in base::names(all_clusterings)){
        PCAs[[base::paste0("module ", a)]] <- plot_PCA_cluster(gtc = all_clusterings[[a]], algo = a, PCA_save_folder = "pca_algo_compare", cols = cols)
      }
    }else{
      PCAs[[base::paste0("module ", algo)]] <- plot_PCA_cluster(gtc = all_clusterings[[algo]], algo = algo, PCA_save_folder = "pca_algo_compare", cols = cols)
    }

  }else{
    PCAs[[base::paste0("module my_clusters")]] <- plot_PCA_cluster(gtc = gtc, algo = "my clusters", PCA_save_folder = "pca_algo_compare", cols = cols)
  }

  hcobject[["integrated_output"]][["PCAs"]] <<- PCAs
}
