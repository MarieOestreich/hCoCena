#' Functional enrichment analysis
#' 
#' Returns a cluster-wise Gene Ontology Database enrichment.
#' @param gene_sets A vector. The names of databases enrichment should be performed for. Choose one or multiple of "Go", "Kegg", "Hallmark", and/or "Reactome".
#'  Available databases depend on supplement files previously set
#' @param top Integer. The number of most strongly enriched terms to return per cluster. Default is 5.
#' @param clusters Either "all" (default) or a vector of clusters as strings. Defines for which clusters to perform the enrichment.
#' @param padj Method to use for multiple testing correction. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".  Default is "BH" (Benjamini-Hochberg).
#' @param qval Upper threshold for the adjusted p-value. Default is 0.05.
#' @export

functional_enrichment <- function(gene_sets = c("Go", "Kegg", "Hallmark", "Reactome"), top = 5,  clusters = c("all"), padj = "BH", qval = 0.05){
  
  # Define general stats
  universe <- base::lapply(1:base::length(hcobject[["layers"]]), 
                           function(x) {
                             return(base::rownames(hcobject[["data"]][[base::paste0("set", x, "_counts")]]))
                           }) %>% base::unlist() %>% base::unique()
  
  cluster_info <- hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]
  
  all_clusters <- base::unique(cluster_info$color)
  all_clusters <- all_clusters[all_clusters != "white"]
  
  if (clusters[1] == "all") {
    clusters <- all_clusters
  }

  # Perform the analysis for each of the selected gene sets
  for (i in gene_sets){
    if(i %in% names(hcobject[["supplementary_data"]])){
      
      top_enr <- list()
      res <- list()
      
      # Perform the enrichment for each of the clusters
      for (c in clusters) {
        
        genes <- dplyr::filter(cluster_info, color == c) %>% dplyr::pull(., "gene_n") %>% base::strsplit(., split = ",") %>% base::unlist(.)
        
        enrich <- clusterProfiler::enricher(gene = genes, TERM2GENE = hcobject[["supplementary_data"]][[stringr::str_to_title(i)]], 
                                            qvalueCutoff = qval, pAdjustMethod = padj, universe = universe)
        
        if (!base::is.null(enrich)) {
          tmp <- enrich@result
          tmp <- dplyr::filter(tmp, qvalue <= qval)
          if (base::nrow(tmp) == 0) {
            next
          } else {
            tmp <- tmp[base::order(tmp$qvalue, decreasing = F), ]
            top_enr[[c]] <- tmp$Description[1:top]
            res[[c]] <- enrich@result
          }
        }
      }
      
      top_enr <- dplyr::bind_rows(top_enr)
      
      if (base::nrow(top_enr) == 1) {
        top_enr <- base::t(top_enr)
      }
      
        
      # Enrichment plot
      plot_order <- ComplexHeatmap::row_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]])
      cluster_order <- all_clusters[plot_order] 
      cluster_order <- cluster_order[cluster_order %in% clusters]
      
      cols <- cluster_order
      names(cols) <- cluster_order
        
      ggplot_df <- base::data.frame(terms   = base::as.vector(base::as.matrix(top_enr)), 
                                    cluster = base::factor(base::rep(base::colnames(top_enr), each = base::nrow(top_enr)), levels = cluster_order), 
                                    val     = base::factor(base::rep(base::colnames(top_enr), each = base::nrow(top_enr)), levels = cluster_order) %>% base::as.numeric())
      ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), ] %>% dplyr::arrange(cluster)
      ggplot_df$terms <- base::factor(ggplot_df$terms, levels = base::unique(ggplot_df$terms))
      
      p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, y = terms, fill = cluster)) + 
        ggplot2::geom_vline(xintercept = base::c(1:base::length(cluster_order)), color = cluster_order, alpha = 0.5, size = 1.5) + 
        ggplot2::geom_point(size = 6, pch = 21, color = "white", stroke = 1.5) + 
        ggplot2::scale_fill_manual(values = cols) + 
        ggplot2::scale_x_discrete(limits = cluster_order) +
        ggplot2::scale_y_discrete(position = "right") + 
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                       axis.ticks.x = ggplot2::element_blank(), 
                       legend.position = "none", 
                       panel.grid.major.x = ggplot2::element_blank(), 
                       panel.grid.major.y = ggplot2::element_line(size = 0.1, color = "lightgrey")) + 
        ggplot2::xlab("") + 
        ggplot2::ylab("") + 
        ggplot2::ggtitle(paste0(stringr::str_to_title(i), " enrichment"))
      
      # Heatmap 
      col_order <- ComplexHeatmap::column_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]])
      m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
      mm <- reshape2::melt(m)
      mm$Var2 <- base::factor(mm$Var2, levels = base::unique(mm$Var2)[col_order])
      
      p2 <- ggplot2::ggplot(mm, ggplot2::aes(x = Var1, y = Var2, fill = value)) + 
        ggplot2::geom_tile(color = "white", height = 1) + 
        ggplot2::scale_fill_gradientn(colors = (grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))))(base::length(base::seq(-2, 2, by = 0.1)))) + 
        ggplot2::scale_x_discrete(limits = cluster_order) +
        ggplot2::scale_y_discrete(position = "right") + 
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
                       legend.position = "none") + 
        ggplot2::xlab("") +
        ggplot2::ylab("")
  
      
      ## Combine the plots and format the output
      
      cp <- egg::ggarrange(p, p2, ncol = 1, heights = c(3, 1))
      Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Enrichment_", stringr::str_to_title(i), "_top_", top, ".pdf"), 
                      width = 8, height = 8)
      print(cp)
      grDevices::dev.off()
      openxlsx::write.xlsx(x = res, 
                           file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Enrichment_", stringr::str_to_title(i), ".xlsx"), 
                           overwrite = T)
      output <- list()
      output[["p"]] <- cp
      output[["result"]] <- top_enr
      output[["enrichment"]] <- res
      hcobject[["integrated_output"]][["enrichments"]][[paste0("top_", i)]] <<- output
      
    } else {
      print(paste0("invalid database: ",i))
    }
  }
}
