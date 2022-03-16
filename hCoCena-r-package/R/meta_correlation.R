#' Numeric Metadata Correaltion
#' 
#' Calculates a matrix where columns are names of the meta categories and rows are clusters. Cells contain the Pearson correlation value between 
#' 		a) the mean expressions of cluster genes in each sample with 
#' 		b) the numeric meta value in each sample.
#' @param set An integer. The number of the datset for which to perform the calculation.
#' @param meta A vector of strings. The names of the numeric annotation column(s) which to correlate to the cluster expression patterns.
#' @param p_val The maximum p-value to determine a correlation as significant. Default is 0.05. Non-significant correlations are shown in grey.
#' @export


meta_correlation_num <- function(set, meta, p_val = 0.05){
  
  meta_data <- dplyr::select(hcobject[["data"]][[base::paste0("set", set, "_anno")]], meta)
  counts <- sample_wise_cluster_expression(set = set)
  cors <- base::lapply(base::colnames(meta_data), function(x){
    meta_vals <- dplyr::pull(meta_data, x) %>% base::as.character() %>% base::as.numeric()
    correlation <- base::lapply(base::rownames(counts), function(y){
      mean_vals <- counts[y,] %>% base::t() %>% base::as.data.frame() %>% dplyr::pull(., y)
      cor_and_pval <- stats::cor.test(mean_vals, meta_vals, method = "pearson") 
      return(base::data.frame(group_x = x, group_y = y, r = (cor_and_pval$estimate %>% base::round(., digits = 2)), p = cor_and_pval$p.value ))
  
    }) %>% rlist::list.rbind() %>% base::as.data.frame()
    return(correlation)
  }) %>% rlist::list.rbind()

  
  cors$color <- base::apply(cors, 1, function(x){
    if(as.double(x["p"]) > p_val){
      return(NA)
    }else{
      return(as.double(x["r"]))
    }
  }) %>% base::as.numeric()

  
  
  g <- ggplot2::ggplot(data = cors, ggplot2::aes(x = group_x, 
                               y = base::factor(group_y, levels = base::rev(base::unique(group_y)[ComplexHeatmap::row_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]])])), 
                               fill = color)) + 
    ggplot2::geom_tile(color = "black")+
    ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 7, name = "BrBG")))(21), limits = c(-1,1))+
    ggplot2::geom_text(ggplot2::aes(group_x, group_y, label = r), color = "black", size = 4)+
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90))+
    ggplot2::scale_y_discrete(position = "right")
  
  graphics::plot(g)
  ggplot2::ggsave(base::paste0("numerical_meta_correlation_", hcobject[["layers_names"]][set], ".pdf"),
                  g, device = cairo_pdf, width = 10, height = 8, path = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]]))
  
}




#' Categorical Metadata Correaltion
#' 
#' Calculates as matrix where columns are the different values of the meta category (e.g. if meta = survival and "survival can have values "yes" or "no", then "yes" and "no" will be the columns) and rows are clusters. 
#' 	Cells contain the Pearson correlation value between 
#' 		a) the mean expressions of cluster genes in each voi group with 
#' 		b) the counts of the meta value (e.g., "yes" or "no") across voi groups.
#' @param set An integer. The number of the datset for which to perform the calculation.
#' @param meta A single string. The name of the categorical annotation column which to correlate to the cluster expression patterns.
#' @param p_val The maximum p-value to determine a correlation as significant. Default is 0.05. Non-significant correlations are shown in grey.
#' @export

meta_correlation_cat <- function(meta, set, p_val = 0.05){
  
  # Correlation of the pattern the metainfo has across groups with the pattern of the clusters across groups 
  # can only be calculated when there are at least 2 groups
  if(base::ncol(hcobject[["layer_specific_outputs"]][[base::paste0("set", set)]][["part2"]][["GFC_all_genes"]]) == 2){
    stop("More than one group is needed to caclulate a correlation. The only group present is:    ",
         base::colnames(hcobject[["layer_specific_outputs"]][[base::paste0("set", set)]][["part2"]][["GFC_all_genes"]])[1])
  }
  
  #  extract anno of given data set:
  anno <- hcobject[["data"]][[base::paste0("set", set, "_anno")]]
  
  # for convenience: data frame of 3 columns - one containing sample IDs, one with corresponding voi-value and one with
  # corresponding value of the meta information of interest
  df <- base::data.frame(sample = base::rownames(anno), 
                   voi = anno[hcobject[["global_settings"]][["voi"]]], 
                   meta = anno[meta])
  base::colnames(df) <- base::c("sample", "voi", "meta")

  
  # create data frame with clusters as rows and samples as columns. The value of a cell [i,j] is the mean expression of
  # all genes in cluster i in sample j:
  cluster_X_sample_counts <- sample_wise_cluster_expression(set = set)
  
  # Group the samples into their voi groups to correspond to the module heatmap. 
  # Grouping of values happens by calculating their mean:
  cluster_X_sample_mean_counts <- base::lapply(base::unique(df$voi), function(x){
    samples <- dplyr::filter(df, voi == x) %>% dplyr::pull(., "sample")
    tmp <- dplyr::select(cluster_X_sample_counts, samples) %>% base::apply(., 1, base::mean)
    tmp <- base::data.frame(v1 = tmp)
    base::colnames(tmp) <- x
    return(tmp)
  }) %>% rlist::list.cbind()
  
#  # Remove control group, if control is not "none":
#  cluster_X_sample_mean_counts <- cluster_X_sample_mean_counts[, base::colnames(cluster_X_sample_mean_counts) %in%
#                                                                 base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])]

  
  # For each group of samples (voi) store the values of the meta info of interest for these samples:
  voi_list <- base::lapply(base::unique(df$voi), function(x){
    tmp <- dplyr::filter(df, voi == x)
    return(tmp["meta"])
  })
  
  base::names(voi_list) <- base::unique(df$voi)
  
  # For each value that the meta information can take (e.g., recovered = YES or NO) count its occurances per sample group:
  out <- base::lapply(base::unique(df$meta), function(x){
    counts <- base::lapply(voi_list, function(y){
      vec <- y[,1]
      return(base::length(vec[vec == x]))
    }) %>% base::unlist()
    out <- base::data.frame(counts = counts)
    base::rownames(out) <- base::names(voi_list)
    base::colnames(out) <- base::c(base::as.character(x))

    return(out)
  }) 
  out <- rlist::list.cbind(out) %>% base::t() %>% base::as.data.frame()
  out <- dplyr::select(out, base::colnames(cluster_X_sample_mean_counts))

  
  # Transform absolute occurrance counts to percentages:
  out <- base::apply(out, 1, function(x){
    if(base::sum(x) == 0){
      base::rep(0, base::length(x))
    }else{
      x/base::sum(x)
    }
  })%>% base::t() %>% base::as.data.frame()

  
  
  # Calculate pearson correlation of those count fractions across groups with cluster mean expressions
  # across groups:
  cors <- base::lapply(1:base::nrow(out), function(x){
    vec1 <- base::as.numeric(out[x,])
    tmp <- base::lapply(1:base::nrow(cluster_X_sample_mean_counts), function(y){
      vec2 <- base::as.numeric(cluster_X_sample_mean_counts[y,])
      cor_res <- stats::cor.test(x = vec1, y = vec2, method = "pearson")
      return(base::data.frame(r = cor_res$estimate, p = cor_res$p.value))
    }) %>% rlist::list.rbind()
    base::rownames(tmp) <- base::rownames(cluster_X_sample_mean_counts)
    return(tmp)
  })
  base::names(cors) <- base::rownames(out)
  
  
  
  # Prepare for plotting:
  vals <- base::lapply(base::names(cors), function(x){
    tmp <- base::cbind(base::data.frame(group_x = base::rep(x, base::nrow(cors[[x]])), group_y = base::rownames(cors[[x]])), cors[[x]])
    return(tmp)
  }) %>% rlist::list.rbind()
  
  vals$color <- base::apply(vals, 1, function(x){
    if(x["p"] > p_val){
      return(NA)
    }else{
      return(x["r"])
    }
  }) %>% base::as.numeric()
  
  vals["r"] <- base::round(vals["r"], digits = 2)
  
  
  
  
  # plot
  g <- ggplot2::ggplot(data = vals, ggplot2::aes(x = group_x, 
                          y = base::factor(group_y, levels = base::rev(base::unique(group_y)[ComplexHeatmap::row_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]])])), 
                          fill = color)) + 
    ggplot2::geom_tile(color = "black")+
    ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 7, name = "BrBG")))(21), limits = c(-1,1))+
    ggplot2::geom_text(ggplot2::aes(group_x, group_y, label = r), color = "black", size = 4)+
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank())+
    ggplot2::scale_y_discrete(position = "right")
  
  graphics::plot(g)
  
}
