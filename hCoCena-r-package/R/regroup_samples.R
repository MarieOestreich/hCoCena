#' Cut Hierarchical Clustering Tree For Sample Regrouping
#' 
#' If the variable of interest (voi) does not go well with the clustering of the heatmaps returned by "run_expression_analysis_2" or as seen in the PCA, 
#' 	new group labels can be assigned to the samples based on the data structure rather than meta information. 
#' 	This allows you to analyse the data by defining unknown subgroups and thus not solely rely on prior knowledge. 
#' 	The regrouping can be performed using hierarchical clustering on all genes, the network genes or the modules.
#' 	Note: This function can be rerun to try out multiple settings for the k-parameter, if the 'save' parameter is set to FALSE. If it is set to TRUE, the original groups will be overwritten.
#' 	Thus, only set 'save' to TRUE once you have decided on your k to cut the tree.
#' @param by This parameter defines based on what data to regroup the samples. It can be one of: "all" (the gene expressions of all genes are used (DEFAULT)), 
#' 	"network" (only the expressions of genes present in the layer's network are used) or "module" (the regrouping is done based on the samples' mean expression patterns across modules).
#' 	Please note that in case of by = "network" and by = "module" the output is a heatmap including a dendrogram, while in the case of by = "all" the output is only the dendrogram. 
#' 	This is due to performance reasons, since heatmaps spanning the entire gene space are often to memory intensive to plot.
#' @param set Defines the layers/datasets of which the samples shall be regrouped. Either "all" (provides a hierarchical clustering for all datasets based on which samples can be regrouped (DEFAULT))
#' 	or as vector of integers e.g. c(1,3,...), only the layers corresponding to the numbers given in the vector are up for regrouping.
#' @param method The method used for the agglomerative clustering, for options refer to the 'hclust()' documentation, DEFAULT is 'complete'.
#' @param k A vector of integers giving for each layer that is meant to be regrouped the number of clusters to retrieve after cutting the dendrogram. 
#' 	Default is a vector of 1s (no cutting, good for first inspection).
#' @param save A Boolean, default is FALSE. 
#' 	Once you have decided on the k-values for cutting the trees and regrouping the samples accordingly, RE-RUN the function, and set this parameter to TRUE. 
#' 	This will export a new annotation file where the (voi) column now contains the new groups and it will also alter you annotation object in hcobject$data accordingly 
#' 	(you will find the original groups in a new column labelled (voi)_old ). 
#' 	The file will be saved in the same directory as the original annotation, bearing the name of the annotation file with the '_regrouped' suffix. 
#' 	THEN RE-RUN CoCena FROM THE 'run_expression_analysis_2' FUNCTION ONWARD (including that function). 
#' 	NOTE: New annotation files are always exported without rownames. 
#' @export


cut_hclust <- function(by = "all", set = "all", method = "complete", k = base::rep(1, base::length(hcobject[["layers"]])), save = F){
  # plots heatmaps AND dendrograms for options "module" and "network", but only dendrogram for option "all". 
  # This is due to performance issues with massive heatmaps like that

  if(!by %in% base::c("all", "network", "module")){
    stop("'", by, "' is not a valid value for the parameter 'by'. Choose one of 'all', 'network' or 'module'")
  }
  
  if(by == "all"){
    if(set == "all"){
      for(l in 1:base::length(hcobject[["layers"]])){
        dat <- hcobject[["data"]][[base::paste0("set", l, "_counts")]] %>% 
        	base::as.matrix()
        
        
        dd <- stats::dist(base::t(dat), method = "euclidean")
        hc <- stats::hclust(dd, method = method)
        cut <- stats::cutree(hc, k = k[l])
        dend <- stats::as.dendrogram(hc)
        
        dend <- dend %>% dendextend::set("branches_k_color", 
                                         value = ggsci::pal_d3(palette = "category20")(20)[1:k[l]], 
                                         k = k[l]) 
        graphics::plot(dend)
        
        if(save){
          new_group <- base::data.frame(sample = base::names(cut), group = cut)
          base::rownames(new_group) <- new_group$sample
          new_group <- new_group[base::rownames(hcobject[["data"]][[base::paste0("set", l, "_anno")]]),]
          anno <- hcobject[["data"]][[base::paste0("set", l, "_anno")]]
          anno[[base::paste0(hcobject[["global_settings"]][["voi"]], "_old")]] <- anno[[hcobject[["global_settings"]][["voi"]]]]
          anno[[hcobject[["global_settings"]][["voi"]]]] <- base::paste0(hcobject[["layers_names"]][l], "_", new_group$group)
          filename <- base::strsplit(hcobject[["layers"]][[base::paste0("set", l)]][2], split = ".txt") %>% 
          	base::unlist()
          if(! hcobject[["global_settings"]][["control"]] == "none"){
            anno <- find_new_ctrl(anno = anno, l = l)
          }
          utils::write.table(x = anno, 
                      file = base::paste0(hcobject[["working_directory"]][["dir_annotation"]], filename, "_regrouped.txt"), 
                      sep = "\t", quote = F, row.names = F, col.names = T)
          hcobject[["data"]][[base::paste0("set", l, "_anno")]] <<- anno
          message("After the regrouping, you must rerun the steps from the main markdown beginning with run_expression_analysis_2(), since the sample groups have now changed.")
        }
      } # end for loop
    }else{
      for(l in set){
        dat <- hcobject[["data"]][[base::paste0("set", l, "_counts")]] %>% 
        	base::as.matrix()
        
        dd <- stats::dist(base::t(dat), method = "euclidean")
        hc <- stats::hclust(dd, method = method)
        
        cut <- stats::cutree(hc, k = k[l])
        dend <- stats::as.dendrogram(hc)
        
        dend <- dend %>% dendextend::set("branches_k_color", 
                                         value = ggsci::pal_d3(palette = "category20")(20)[1:k[l]], 
                                         k = k[l]) 
        graphics::plot(dend)
        if(save){
          new_group <- base::data.frame(sample = base::names(cut), group = cut)
          base::rownames(new_group) <- new_group$sample
          new_group <- new_group[base::rownames(hcobject[["data"]][[base::paste0("set", l, "_anno")]]),]
          anno <- hcobject[["data"]][[base::paste0("set", l, "_anno")]]
          anno[[base::paste0(hcobject[["global_settings"]][["voi"]], "_old")]] <- anno[[hcobject[["global_settings"]][["voi"]]]]
          anno[[hcobject[["global_settings"]][["voi"]]]] <- base::paste0(hcobject[["layers_names"]][l], "_", new_group$group)
          filename <- base::strsplit(hcobject[["layers"]][[base::paste0("set", l)]][2], split = ".txt") %>% base::unlist()
          if(! hcobject[["global_settings"]][["control"]] == "none"){
            anno <- find_new_ctrl(anno = anno, l = l)
          }
          utils::write.table(x = anno, 
                      file = paste0(hcobject[["working_directory"]][["dir_annotation"]], filename, "_regrouped.txt"), 
                      sep = "\t", quote = F, row.names = F, col.names = T)
          hcobject[["data"]][[base::paste0("set", l, "_anno")]] <<- anno
          message("After the regrouping, you must rerun the steps from the main markdown beginning with run_expression_analysis_2(), since the sample groups have now changed.")
        }
      } # end for loop
    } # end if-else
  } # end by == "all"
  
  if(by == "network"){
    
    if(set == "all"){
      for(l in 1:base::length(hcobject[["layers"]])){
        
        if(! base::paste0("set", l) %in% base::names(hcobject[["layer_specific_outputs"]])){
          stop("You have to have created layer-specific networks to perform this step.")
        }
        if(! "part2" %in% base::names(hcobject[["layer_specific_outputs"]][[base::paste0("set", l)]]) ){
          stop("You have to have created layer-specific networks to perform this step.")
        }
        if(! "heatmap_out" %in% base::names(hcobject[["layer_specific_outputs"]][[base::paste0("set", l)]][["part2"]]) ){
          stop("You have to have created layer-specific networks to perform this step.")
        }
        if(! "filt_cutoff_graph" %in% base::names(hcobject[["layer_specific_outputs"]][[base::paste0("set", l)]][["part2"]][["heatmap_out"]]) ){
          stop("You have to have created layer-specific networks to perform this step.")
        }
 
        
        dat <- hcobject[["data"]][[base::paste0("set", l, "_counts")]] %>% base::as.matrix()
        genes <- igraph::V(hcobject[["layer_specific_outputs"]][[base::paste0("set", l)]][["part2"]][["heatmap_out"]][["filt_cutoff_graph"]])$name 
        dat <- dat[genes,]
        p <- pheatmap::pheatmap(mat = dat,
                                color = base::rev(RColorBrewer::brewer.pal(11, "RdBu")),
                                scale = "row",
                                cluster_rows =  T,
                                cluster_cols =  T,
                                fontsize = 8,
                                show_rownames = F, 
                                show_colnames = T, 
                                clustering_distance_cols = "euclidean", 
                                clustering_method = method, 
                                cutree_cols = k[l])
        
        
        p
        if(save){
          cut <- stats::cutree(p$tree_col, k=k[l])
          new_group <- base::data.frame(sample = base::names(cut), group = cut)
          base::rownames(new_group) <- new_group$sample
          new_group <- new_group[base::rownames(hcobject[["data"]][[base::paste0("set", l, "_anno")]]),]
          anno <- hcobject[["data"]][[base::paste0("set", l, "_anno")]]
          anno[[base::paste0(hcobject[["global_settings"]][["voi"]], "_old")]] <- anno[[hcobject[["global_settings"]][["voi"]]]]
          anno[[hcobject[["global_settings"]][["voi"]]]] <- base::paste0(hcobject[["layers_names"]][l], "_", new_group$group)
          filename <- base::strsplit(hcobject[["layers"]][[base::paste0("set", l)]][2], split = ".txt") %>% 
          	base::unlist()
          if(! hcobject[["global_settings"]][["control"]] == "none"){
            anno <- find_new_ctrl(anno = anno, l = l)
          }
          utils::write.table(x = anno, 
                      file = base::paste0(hcobject[["working_directory"]][["dir_annotation"]], filename, "_regrouped.txt"), 
                      sep = "\t", quote = F, row.names = F, col.names = T)
          hcobject[["data"]][[base::paste0("set", l, "_anno")]] <<- anno
          message("After the regrouping, you must rerun the steps from the main markdown beginning with run_expression_analysis_2(), since the sample groups have now changed.")
        }
      } # end for loop
    }else{
      for(l in set){
        if(! base::paste0("set", l) %in% base::names(hcobject[["layer_specific_outputs"]])){
          stop("You have to have created layer-specific networks to perform this step.")
        }
        if(! "part2" %in% base::names(hcobject[["layer_specific_outputs"]][[base::paste0("set", l)]]) ){
          stop("You have to have created layer-specific networks to perform this step.")
        }
        if(! "heatmap_out" %in% base::names(hcobject[["layer_specific_outputs"]][[base::paste0("set", l)]][["part2"]]) ){
          stop("You have to have created layer-specific networks to perform this step.")
        }
        if(! "filt_cutoff_graph" %in% base::names(hcobject[["layer_specific_outputs"]][[base::paste0("set", l)]][["part2"]][["heatmap_out"]]) ){
          stop("You have to have created layer-specific networks to perform this step.")
        }
        dat <- hcobject[["data"]][[base::paste0("set", l, "_counts")]] %>% base::as.matrix()
        genes <- igraph::V(hcobject[["layer_specific_outputs"]][[base::paste0("set", l)]][["part2"]][["heatmap_out"]][["filt_cutoff_graph"]])$name 
        dat <- dat[genes,]
        p <- pheatmap::pheatmap(mat = dat,
                                color=base::rev(RColorBrewer::brewer.pal(11, "RdBu")),
                                scale="row",
                                cluster_rows=T,
                                cluster_cols=T,
                                fontsize = 8,
                                show_rownames = F, 
                                show_colnames = T, 
                                clustering_distance_cols="euclidean", 
                                clustering_method=method, 
                                cutree_cols = k[l])
        p
        if(save){
          cut <- stats::cutree(p$tree_col, k=k[l])
          new_group <- base::data.frame(sample = base::names(cut), group = cut)
          base::rownames(new_group) <- new_group$sample
          new_group <- new_group[base::rownames(hcobject[["data"]][[base::paste0("set", l, "_anno")]]),]
          anno <- hcobject[["data"]][[base::paste0("set", l, "_anno")]]
          anno[[base::paste0(hcobject[["global_settings"]][["voi"]], "_old")]] <- anno[[hcobject[["global_settings"]][["voi"]]]]
          anno[[hcobject[["global_settings"]][["voi"]]]] <- base::paste0(hcobject[["layers_names"]][l], "_", new_group$group)
          filename <- base::strsplit(hcobject[["layers"]][[base::paste0("set", l)]][2], split = ".txt") %>% base::unlist()
          if(! hcobject[["global_settings"]][["control"]] == "none"){
            anno <- find_new_ctrl(anno = anno, l = l)
          }
          utils::write.table(x = anno, 
                      file = base::paste0(hcobject[["working_directory"]][["dir_annotation"]], filename, "_regrouped.txt"), 
                      sep = "\t", quote = F, row.names = F, col.names = T)
          hcobject[["data"]][[base::paste0("set", l, "_anno")]] <<- anno
          message("After the regrouping, you must rerun the steps from the main markdown beginning with run_expression_analysis_2(), since the sample groups have now changed.")
        }
      } # end for loop
    } # end if-else
  } # end by == "network"
  
  
  if(by == "module"){
    if(base::is.null(hcobject[["integrated_output"]])){
      stop("You have to have clustered the network to perform this step.")
    }
    if(! "cluster_calc" %in% base::names(hcobject[["integrated_output"]])){
      stop("You have to have clustered the network to perform this step.")
    }
    if(! "cluster_information" %in% base::names(hcobject[["integrated_output"]][["cluster_calc"]])){
      stop("You have to have clustered the network to perform this step.")
    }
    if(set == "all"){
      for(l in 1:base::length(hcobject[["layers"]])){
        
        
        gtc <- GeneToCluster()
        counts <- hcobject[["data"]][[base::paste0("set", l, "_counts")]]
        FCs <- NULL
        for(c in base::unique(gtc$color)){
          if(!c == "white"){
            genes <- gtc[gtc$color == c,] %>% dplyr::pull(., "gene")
            tmp_counts <- counts[base::rownames(counts) %in% genes, ]
            tmp_fc <- base::apply(tmp_counts, 1, function(x){
              vec <- gtools::foldchange(x, base::mean(x))
              vec_no_inf <- vec[!vec == Inf & !vec == -Inf]
              vec[vec == Inf] <- base::max(vec_no_inf)
              vec[vec == -Inf] <- base::min(vec_no_inf)
              return(vec)
            }) %>% base::t() %>% base::as.data.frame()
            tmp_fc <- base::apply(tmp_fc, 2, base::mean) %>% base::t() %>% base::as.data.frame()
            base::rownames(tmp_fc) <- c
            FCs <- base::rbind(FCs, tmp_fc)
          }
        }
        base::colnames(FCs) <- base::colnames(counts)
        FCs[FCs > hcobject[["global_settings"]][["range_GFC"]]] <- hcobject[["global_settings"]][["range_GFC"]]
        FCs[FCs < -hcobject[["global_settings"]][["range_GFC"]]] <- -hcobject[["global_settings"]][["range_GFC"]]
        # saving FCs of sample X cluster:
        hcobject[["satellite_outputs"]][["FC_samples_X_clusters"]][[base::paste0("layer_", l)]] <<- FCs
        
        p <- pheatmap::pheatmap(mat = FCs,
                                color=grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(base::length(base::seq(-2, 2, by = .1))),
                                scale="column",
                                cluster_rows=F,
                                cluster_cols=T,
                                fontsize = 8,
                                show_rownames = T, 
                                show_colnames = T, 
                                clustering_distance_cols="euclidean", 
                                clustering_method=method,
                                breaks = base::seq(-2, 2, by = 0.1), 
                                cutree_cols = k[l])
        p
        if(save){
          cut <- stats::cutree(p$tree_col, k=k[l])
          new_group <- base::data.frame(sample = base::names(cut), group = cut)
          base::rownames(new_group) <- new_group$sample
          new_group <- new_group[base::rownames(hcobject[["data"]][[base::paste0("set", l, "_anno")]]),]
          anno <- hcobject[["data"]][[base::paste0("set", l, "_anno")]]
          anno[[base::paste0(hcobject[["global_settings"]][["voi"]], "_old")]] <- anno[[hcobject[["global_settings"]][["voi"]]]]
          anno[[hcobject[["global_settings"]][["voi"]]]] <- base::paste0(hcobject[["layers_names"]][l], "_", new_group$group)
          filename <- base::strsplit(hcobject[["layers"]][[base::paste0("set", l)]][2], split = ".txt") %>% base::unlist()
          if(! hcobject[["global_settings"]][["control"]] == "none"){
            anno <- find_new_ctrl(anno = anno, l = l)
          }
          utils::write.table(x = anno, 
                      file = base::paste0(hcobject[["working_directory"]][["dir_annotation"]], filename, "_regrouped.txt"), 
                      sep = "\t", quote = F, row.names = F, col.names = T)
          hcobject[["data"]][[base::paste0("set", l, "_anno")]] <<- anno
          message("After the regrouping, you must rerun the steps from the main markdown beginning with run_expression_analysis_2(), since the sample groups have now changed.")
        }
      } # end for loop
    }else{
      
      for(l in set){
        
        gtc <- GeneToCluster()
        counts <- hcobject[["data"]][[base::paste0("set", l, "_counts")]]
        FCs <- NULL
        for(c in base::unique(gtc$color)){
          if(!c == "white"){
            genes <- gtc[gtc$color == c,] %>% dplyr::pull(., "gene")
            tmp_counts <- counts[base::rownames(counts) %in% genes, ]
            tmp_fc <- base::apply(tmp_counts, 1, function(x){
              vec <- gtools::foldchange(x, base::mean(x))
              vec_no_inf <- vec[!vec == Inf & !vec == -Inf]
              vec[vec == Inf] <- base::max(vec_no_inf)
              vec[vec == -Inf] <- base::min(vec_no_inf)
              return(vec)
            }) %>% base::t() %>% base::as.data.frame()
            tmp_fc <- base::apply(tmp_fc, 2, base::mean) %>% base::t() %>% base::as.data.frame()
            base::rownames(tmp_fc) <- c
            FCs <- base::rbind(FCs, tmp_fc)
          }
        }
        base::colnames(FCs) <- base::colnames(counts)
        FCs[FCs > hcobject[["global_settings"]][["range_GFC"]]] <- hcobject[["global_settings"]][["range_GFC"]]
        FCs[FCs < -hcobject[["global_settings"]][["range_GFC"]]] <- -hcobject[["global_settings"]][["range_GFC"]]
        # saving FCs of sample X cluster:
        hcobject[["satellite_outputs"]][["FC_samples_X_clusters"]][[base::paste0("layer_", l)]] <<- FCs
        
        p <- pheatmap::pheatmap(mat = FCs,
                                color=grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(base::length(base::seq(-2, 2, by = .1))),
                                scale="column",
                                cluster_rows=F,
                                cluster_cols=T,
                                fontsize = 8,
                                show_rownames = T, 
                                show_colnames = T, 
                                clustering_distance_cols="euclidean", 
                                clustering_method=method,
                                breaks = base::seq(-2, 2, by = 0.1), 
                                cutree_cols = k[l])
        p
        if(save){
          cut <- stats::cutree(p$tree_col, k=k[l])
          new_group <- base::data.frame(sample = base::names(cut), group = cut)
          base::rownames(new_group) <- new_group$sample
          new_group <- new_group[base::rownames(hcobject[["data"]][[base::paste0("set", l, "_anno")]]),]
          anno <- hcobject[["data"]][[base::paste0("set", l, "_anno")]]
          anno[[base::paste0(hcobject[["global_settings"]][["voi"]], "_old")]] <- anno[[hcobject[["global_settings"]][["voi"]]]]
          anno[[hcobject[["global_settings"]][["voi"]]]] <- base::paste0(hcobject[["layers_names"]][l], "_", new_group$group)
          filename <- base::strsplit(hcobject[["layers"]][[base::paste0("set", l)]][2], split = ".txt") %>% base::unlist()
          if(! hcobject[["global_settings"]][["control"]] == "none"){
            anno <- find_new_ctrl(anno = anno, l = l)
          }
          utils::write.table(x = anno, 
                      file = base::paste0(hcobject[["working_directory"]][["dir_annotation"]], filename, "_regrouped.txt"), 
                      sep = "\t", quote = F, row.names = F, col.names = T)
          hcobject[["data"]][[base::paste0("set", l, "_anno")]] <<- anno
          message("After the regrouping, you must rerun the steps from the main markdown beginning with run_expression_analysis_2(), since the sample groups have now changed.")
        }
      } # end for loop
    } # end if-else
  } # end by == "module"
} 