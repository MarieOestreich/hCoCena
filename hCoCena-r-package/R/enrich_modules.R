#' Functional enrichment analysis
#' 
#' Returns a cluster-wise Gene Ontology Database enrichment.
#' @param gene_sets A vector. The names of databases enrichment should be performed for. Choose one or multiple of "GO", "KEGG", "HALLMARK", and/or "Reactome".
#'  Available databases depend on supplement files previously set
#' @param top Integer. The number of most strongly enriched terms to return per cluster.
#' @param clusters Either "all" (default) or a vector of clusters as strings. Defines for which clusters to perform the enrichment.
#' @param padj Method to use for multiple testing correction. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".  Default is "BH" (Benjamini-Hochberg).
#' @param qval Upper threshold for the adjusted p-value. Default is 0.05.
#' @param ordered A boolean. TRUE if enriched terms should be ordered by clusters (default) or FALSE if enriched terms should be ordered alphabetically.
#' @export

enrich_modules <- function(gene_sets, top,  clusters = c("all"), padj = "BH", qval = 0.05, ordered = T){

cluster_info <- hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]
universe <- base::lapply(1:base::length(hcobject[["layers"]]), 
                         function(x) {
                           return(base::rownames(hcobject[["data"]][[base::paste0("set", 
                                                                                  x, "_counts")]]))
                         }) %>% base::unlist() %>% base::unique()
if (clusters[1] == "all") {
  clusters <- base::unique(cluster_info$color)
  clusters <- clusters[!clusters == "white"]
}
top_GO <- list()
res <- list()
if ("Hallmark" %in% gene_sets){
  for (c in clusters) {
    genes <- dplyr::filter(cluster_info, color == c) %>% 
      dplyr::pull(., "gene_n") %>% base::strsplit(., split = ",") %>% 
      base::unlist(.)
    enrich <- clusterProfiler::enricher(genes, TERM2GENE = hcobject[["supplementary_data"]][["hallmark"]], 
                                        qvalueCutoff = qval, pAdjustMethod = padj, universe = universe)
    if (!base::is.null(enrich)) {
      tmp <- enrich@result
      tmp <- dplyr::filter(tmp, qvalue <= qval)
      if (base::nrow(tmp) == 0) {
        next
      }
      else {
        tmp <- tmp[base::order(tmp$qvalue, decreasing = F), 
        ]
        if (base::nrow(tmp) < top) {
          diff <- top - base::nrow(tmp)
          top_GO[[c]] <- c(tmp$Description, rep(NA, 
                                                diff))
        }
        else {
          top_GO[[c]] <- tmp$Description[1:top]
        }
        top_GO[[c]] <- tmp$Description[1:top]
        res[[c]] <- enrich@result
      }
    }
  }
  top_GO <- dplyr::bind_rows(top_GO)
  if (base::nrow(top_GO) == 1) {
    top_GO <- base::t(top_GO)
  }
  if(ordered){
    
    plot_order <- ComplexHeatmap::row_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]])
    ggplot_df <- base::data.frame(terms = base::as.vector(base::as.matrix(top_GO)), 
                                  cluster = base::factor(base::rep(base::colnames(top_GO), each = base::nrow(top_GO)), levels = rev(clusters[plot_order])), 
                                  val = base::factor(base::rep(base::colnames(top_GO), 
                                                               each = base::nrow(top_GO)), levels = rev(clusters[plot_order])) %>% base::as.numeric())
    ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), 
    ] %>% dplyr::arrange(desc(cluster))
    ggplot_df$terms <- base::factor(ggplot_df$terms, levels = unique(ggplot_df$terms))
    p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, y = terms)) + 
      ggplot2::geom_vline(xintercept = base::c(1:base::length(base::levels(ggplot_df$cluster))), color = grDevices::adjustcolor(base::as.character(base::levels(ggplot_df$cluster)), alpha.f = 0.5), size = 1.5) + 
      ggplot2::geom_point(ggplot2::aes(fill = cluster), color = "white", stroke = 1.5, size = 6, pch = 21) + 
      ggplot2::scale_fill_manual(values = base::as.character(base::levels(ggplot_df$cluster)[base::levels(ggplot_df$cluster) %in% base::unique(ggplot_df$cluster)])) + 
      ggplot2::scale_color_manual(values = grDevices::adjustcolor(base::as.character(base::levels(ggplot_df$cluster)))) + 
      ggplot2::scale_x_discrete(limits = base::as.character(base::levels(ggplot_df$cluster))) + 
      ggplot2::scale_y_discrete(limits = rev) + ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
                     panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(size = 0.1, 
                                                                                                               color = "lightgrey")) + ggplot2::xlab("") + ggplot2::ylab("") + 
      ggplot2::ggtitle("HALLMARK enrichment")
    m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
    mm <- reshape2::melt(m)
    p2 <- ggplot2::ggplot(mm, ggplot2::aes(Var1, Var2, fill = value)) + 
      ggplot2::geom_tile(color = "white", height = 1) + ggplot2::scale_fill_gradientn(colors = (grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, 
                                                                                                                                                               name = "RdBu"))))(base::length(base::seq(-2, 2, by = 0.1)))) + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::theme_light() + 
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      ggplot2::scale_x_discrete(limits = base::as.character(base::levels(ggplot_df$cluster)))
  }else{      
    ggplot_df <- base::data.frame(terms = base::as.vector(base::as.matrix(top_GO)), 
                                  cluster = base::rep(base::colnames(top_GO), each = base::nrow(top_GO)), 
                                  val = base::factor(base::rep(base::colnames(top_GO), 
                                                               each = base::nrow(top_GO))) %>% base::as.numeric())
    ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), 
    ]
    
    
    
    p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, 
                                                        y = terms)) + ggplot2::geom_vline(xintercept = base::c(1:base::length(base::unique(ggplot_df$cluster))), 
                                                                                          color = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                   "cluster"])), alpha.f = 0.5), size = 1.5) + ggplot2::geom_point(ggplot2::aes(fill = cluster), color = "white", stroke = 1.5, size = 6, pch = 21) + ggplot2::scale_fill_manual(values = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                    "cluster"]))) + ggplot2::scale_color_manual(values = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "cluster"])))) + ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "cluster"]))) + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(size = 0.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       color = "lightgrey")) + ggplot2::xlab("") + ggplot2::ylab("") + 
      ggplot2::ggtitle("HALLMARK enrichment")
    m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
    mm <- reshape2::melt(m)
    p2 <- ggplot2::ggplot(mm, ggplot2::aes(Var1, Var2, fill = value)) + 
      ggplot2::geom_tile(color = "white", height = 1) + ggplot2::scale_fill_gradientn(colors = (grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, 
                                                                                                                                                               name = "RdBu"))))(base::length(base::seq(-2, 2, by = 0.1)))) + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::theme_light() + 
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                   "cluster"])))
  }
  
  cp <- egg::ggarrange(p, p2, ncol = 1, heights = c(3, 1))
  Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], 
                                      hcobject[["global_settings"]][["save_folder"]], "/DJplot_HALLMARK_top_", 
                                      top, ".pdf"), width = 12, height = 12)
  print(cp)
  grDevices::dev.off()
  openxlsx::write.xlsx(x = res, file = base::paste0(hcobject[["working_directory"]][["dir_output"]], 
                                                    hcobject[["global_settings"]][["save_folder"]], "/hallmark_enrichment.xlsx"), 
                       overwrite = T)
  output <- list()
  output[["p"]] <- cp
  output[["result"]] <- top_GO
  output[["enrichment"]] <- res
  hcobject[["integrated_output"]][["enrichments"]][["top_HALL"]] <<- output
  
}
if ("Go" %in% gene_sets){
  top_GO <- list()
  res <- list()
  for (c in clusters) {
    genes <- dplyr::filter(cluster_info, color == c) %>% 
      dplyr::pull(., "gene_n") %>% base::strsplit(., split = ",") %>% 
      base::unlist(.)
    enrich <- clusterProfiler::enricher(genes, TERM2GENE = hcobject[["supplementary_data"]][["go"]], 
                                        qvalueCutoff = qval, pAdjustMethod = padj, universe = universe)
    if (!base::is.null(enrich)) {
      tmp <- enrich@result
      tmp <- dplyr::filter(tmp, qvalue <= qval)
      if (base::nrow(tmp) == 0) {
        next
      }
      else {
        tmp <- tmp[base::order(tmp$qvalue, decreasing = F), 
        ]
        if (base::nrow(tmp) < top) {
          diff <- top - base::nrow(tmp)
          top_GO[[c]] <- c(tmp$Description, rep(NA, 
                                                diff))
        }
        else {
          top_GO[[c]] <- tmp$Description[1:top]
        }
        top_GO[[c]] <- tmp$Description[1:top]
        res[[c]] <- enrich@result
      }
    }
  }
  top_GO <- dplyr::bind_rows(top_GO)
  if (base::nrow(top_GO) == 1) {
    top_GO <- base::t(top_GO)
  }
  if(ordered){
    
    plot_order <- ComplexHeatmap::row_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]])
    ggplot_df <- base::data.frame(terms = base::as.vector(base::as.matrix(top_GO)), 
                                  cluster = base::factor(base::rep(base::colnames(top_GO), each = base::nrow(top_GO)), levels = rev(clusters[plot_order][clusters[plot_order] %in% colnames(top_GO)])), 
                                  val = base::factor(base::rep(base::colnames(top_GO), 
                                                               each = base::nrow(top_GO)), levels = rev(clusters[plot_order][clusters[plot_order] %in% colnames(top_GO)])) %>% base::as.numeric())
    ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), 
    ]%>% dplyr::arrange(desc(cluster))
    ggplot_df$terms <- base::factor(ggplot_df$terms, levels = unique(ggplot_df$terms))
    
    p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, 
                                                        y = terms)) + ggplot2::geom_vline(xintercept = base::c(1:base::length(base::unique(ggplot_df$cluster))), 
                                                                                          color = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                   "cluster"])), alpha.f = 0.5), size = 1.5) + ggplot2::geom_point(ggplot2::aes(fill = cluster), color = "white", stroke = 1.5, size = 6, pch = 21) + ggplot2::scale_fill_manual(values = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                    "cluster"]))) + ggplot2::scale_color_manual(values = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "cluster"])))) + ggplot2::scale_x_discrete(limits = base::as.character(base::levels(ggplot_df$cluster))) + ggplot2::scale_y_discrete(limits = rev) + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(size = 0.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              color = "lightgrey")) + ggplot2::xlab("") + ggplot2::ylab("") + 
      ggplot2::ggtitle("GO enrichment")
    m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
    mm <- reshape2::melt(m)
    p2 <- ggplot2::ggplot(mm, ggplot2::aes(Var1, Var2, fill = value)) + 
      ggplot2::geom_tile(color = "white", height = 1) + ggplot2::scale_fill_gradientn(colors = (grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, 
                                                                                                                                                               name = "RdBu"))))(base::length(base::seq(-2, 2, by = 0.1)))) + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::theme_light() + 
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      ggplot2::scale_x_discrete(limits = base::as.character(base::levels(ggplot_df$cluster)))
  }else{      
    ggplot_df <- base::data.frame(terms = base::as.vector(base::as.matrix(top_GO)), 
                                  cluster = base::rep(base::colnames(top_GO), each = base::nrow(top_GO)), 
                                  val = base::factor(base::rep(base::colnames(top_GO), 
                                                               each = base::nrow(top_GO))) %>% base::as.numeric())
    ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), 
    ]
    
    
    
    p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, 
                                                        y = terms)) + ggplot2::geom_vline(xintercept = base::c(1:base::length(base::unique(ggplot_df$cluster))), 
                                                                                          color = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                   "cluster"])), alpha.f = 0.5), size = 1.5) + ggplot2::geom_point(ggplot2::aes(fill = cluster), color = "white", stroke = 1.5, size = 6, pch = 21) + ggplot2::scale_fill_manual(values = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                    "cluster"]))) + ggplot2::scale_color_manual(values = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "cluster"])))) + ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "cluster"]))) + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(size = 0.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       color = "lightgrey")) + ggplot2::xlab("") + ggplot2::ylab("") + 
      ggplot2::ggtitle("GO enrichment")
    m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
    mm <- reshape2::melt(m)
    p2 <- ggplot2::ggplot(mm, ggplot2::aes(Var1, Var2, fill = value)) + 
      ggplot2::geom_tile(color = "white", height = 1) + ggplot2::scale_fill_gradientn(colors = (grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, 
                                                                                                                                                               name = "RdBu"))))(base::length(base::seq(-2, 2, by = 0.1)))) + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::theme_light() + 
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                   "cluster"])))
  }
  cp <- egg::ggarrange(p, p2, ncol = 1, heights = c(3, 1))
  Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], 
                                      hcobject[["global_settings"]][["save_folder"]], "/DJplot_GO_top_", 
                                      top, ".pdf"), width = 12, height = 12)
  print(cp)
  grDevices::dev.off()
  openxlsx::write.xlsx(x = res, file = base::paste0(hcobject[["working_directory"]][["dir_output"]], 
                                                    hcobject[["global_settings"]][["save_folder"]], "/go_enrichment.xlsx"), 
                       overwrite = T)
  output <- list()
  output[["p"]] <- cp
  output[["result"]] <- top_GO
  output[["enrichment"]] <- res
  hcobject[["integrated_output"]][["enrichments"]][["top_GO"]] <<- output
  
}
if ("Kegg" %in% gene_sets){
  top_GO <- list()
  res <- list()
  for (c in clusters) {
    genes <- dplyr::filter(cluster_info, color == c) %>% 
      dplyr::pull(., "gene_n") %>% base::strsplit(., split = ",") %>% 
      base::unlist(.)
    enrich <- clusterProfiler::enricher(genes, TERM2GENE = hcobject[["supplementary_data"]][["Kegg"]], 
                                        qvalueCutoff = qval, pAdjustMethod = padj, universe = universe)
    if (!base::is.null(enrich)) {
      tmp <- enrich@result
      tmp <- dplyr::filter(tmp, qvalue <= qval)
      if (base::nrow(tmp) == 0) {
        next
      }
      else {
        tmp <- tmp[base::order(tmp$qvalue, decreasing = F), 
        ]
        if (base::nrow(tmp) < top) {
          diff <- top - base::nrow(tmp)
          top_GO[[c]] <- c(tmp$Description, rep(NA, 
                                                diff))
        }
        else {
          top_GO[[c]] <- tmp$Description[1:top]
        }
        top_GO[[c]] <- tmp$Description[1:top]
        res[[c]] <- enrich@result
      }
    }
  }
  top_GO <- dplyr::bind_rows(top_GO)
  if (base::nrow(top_GO) == 1) {
    top_GO <- base::t(top_GO)
  }
  if(ordered){
    
    plot_order <- ComplexHeatmap::row_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]])
    ggplot_df <- base::data.frame(terms = base::as.vector(base::as.matrix(top_GO)), 
                                  cluster = base::factor(base::rep(base::colnames(top_GO), each = base::nrow(top_GO)), levels = rev(clusters[plot_order][clusters[plot_order] %in% colnames(top_GO)])), 
                                  val = base::factor(base::rep(base::colnames(top_GO), 
                                                               each = base::nrow(top_GO)), levels = rev(clusters[plot_order][clusters[plot_order] %in% colnames(top_GO)])) %>% base::as.numeric())
    ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), 
    ] %>% dplyr::arrange(desc(cluster))
    ggplot_df$terms <- base::factor(ggplot_df$terms, levels = unique(ggplot_df$terms))
    
    p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, 
                                                        y = terms)) + ggplot2::geom_vline(xintercept = base::c(1:base::length(base::unique(ggplot_df$cluster))), 
                                                                                          color = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                   "cluster"])), alpha.f = 0.5), size = 1.5) + ggplot2::geom_point(ggplot2::aes(fill = cluster), color = "white", stroke = 1.5, size = 6, pch = 21) + ggplot2::scale_fill_manual(values = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                    "cluster"]))) + ggplot2::scale_color_manual(values = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "cluster"])))) + ggplot2::scale_x_discrete(limits = base::as.character(base::levels(ggplot_df$cluster))) + ggplot2::scale_y_discrete(limits = rev) + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(size = 0.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              color = "lightgrey")) + ggplot2::xlab("") + ggplot2::ylab("") + 
      ggplot2::ggtitle("KEGG enrichment")
    m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
    mm <- reshape2::melt(m)
    p2 <- ggplot2::ggplot(mm, ggplot2::aes(Var1, Var2, fill = value)) + 
      ggplot2::geom_tile(color = "white", height = 1) + ggplot2::scale_fill_gradientn(colors = (grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, 
                                                                                                                                                               name = "RdBu"))))(base::length(base::seq(-2, 2, by = 0.1)))) + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::theme_light() + 
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      ggplot2::scale_x_discrete(limits = base::as.character(base::levels(ggplot_df$cluster)))
  }else{      
    ggplot_df <- base::data.frame(terms = base::as.vector(base::as.matrix(top_GO)), 
                                  cluster = base::rep(base::colnames(top_GO), each = base::nrow(top_GO)), 
                                  val = base::factor(base::rep(base::colnames(top_GO), 
                                                               each = base::nrow(top_GO))) %>% base::as.numeric())
    ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), 
    ]
    
    
    
    p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, 
                                                        y = terms)) + ggplot2::geom_vline(xintercept = base::c(1:base::length(base::unique(ggplot_df$cluster))), 
                                                                                          color = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                   "cluster"])), alpha.f = 0.5), size = 1.5) + ggplot2::geom_point(ggplot2::aes(fill = cluster), color = "white", stroke = 1.5, size = 6, pch = 21) + ggplot2::scale_fill_manual(values = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                    "cluster"]))) + ggplot2::scale_color_manual(values = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "cluster"])))) + ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "cluster"]))) + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(size = 0.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       color = "lightgrey")) + ggplot2::xlab("") + ggplot2::ylab("") + 
      ggplot2::ggtitle("KEGG enrichment")
    m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
    mm <- reshape2::melt(m)
    p2 <- ggplot2::ggplot(mm, ggplot2::aes(Var1, Var2, fill = value)) + 
      ggplot2::geom_tile(color = "white", height = 1) + ggplot2::scale_fill_gradientn(colors = (grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, 
                                                                                                                                                               name = "RdBu"))))(base::length(base::seq(-2, 2, by = 0.1)))) + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::theme_light() + 
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                   "cluster"])))
  }
  cp <- egg::ggarrange(p, p2, ncol = 1, heights = c(3, 1))
  Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], 
                                      hcobject[["global_settings"]][["save_folder"]], "/DJplot_KEGG_top_", 
                                      top, ".pdf"), width = 12, height = 12)
  print(cp)
  grDevices::dev.off()
  openxlsx::write.xlsx(x = res, file = base::paste0(hcobject[["working_directory"]][["dir_output"]], 
                                                    hcobject[["global_settings"]][["save_folder"]], "/kegg_enrichment.xlsx"), 
                       overwrite = T)
  output <- list()
  output[["p"]] <- cp
  output[["result"]] <- top_GO
  output[["enrichment"]] <- res
  hcobject[["integrated_output"]][["enrichments"]][["top_KEGG"]] <<- output
  
}
if ("Reactome" %in% gene_sets){
  top_GO <- list()
  res <- list()
  for (c in clusters) {
    genes <- dplyr::filter(cluster_info, color == c) %>% 
      dplyr::pull(., "gene_n") %>% base::strsplit(., split = ",") %>% 
      base::unlist(.)
    enrich <- clusterProfiler::enricher(genes, TERM2GENE = hcobject[["supplementary_data"]][["Reactome"]], 
                                        qvalueCutoff = qval, pAdjustMethod = padj, universe = universe)
    if (!base::is.null(enrich)) {
      tmp <- enrich@result
      tmp <- dplyr::filter(tmp, qvalue <= qval)
      if (base::nrow(tmp) == 0) {
        next
      }
      else {
        tmp <- tmp[base::order(tmp$qvalue, decreasing = F), 
        ]
        if (base::nrow(tmp) < top) {
          diff <- top - base::nrow(tmp)
          top_GO[[c]] <- c(tmp$Description, rep(NA, 
                                                diff))
        }
        else {
          top_GO[[c]] <- tmp$Description[1:top]
        }
        top_GO[[c]] <- tmp$Description[1:top]
        res[[c]] <- enrich@result
      }
    }
  }
  top_GO <- dplyr::bind_rows(top_GO)
  if (base::nrow(top_GO) == 1) {
    top_GO <- base::t(top_GO)
  }
  if(ordered){
    
    plot_order <- ComplexHeatmap::row_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]])
    ggplot_df <- base::data.frame(terms = base::as.vector(base::as.matrix(top_GO)), 
                                  cluster = base::factor(base::rep(base::colnames(top_GO), each = base::nrow(top_GO)), levels = rev(clusters[plot_order][clusters[plot_order] %in% colnames(top_GO)])), 
                                  val = base::factor(base::rep(base::colnames(top_GO), 
                                                               each = base::nrow(top_GO)), levels = rev(clusters[plot_order][clusters[plot_order] %in% colnames(top_GO)])) %>% base::as.numeric())
    ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), 
    ] %>% dplyr::arrange(desc(cluster))
    ggplot_df$terms <- base::factor(ggplot_df$terms, levels = unique(ggplot_df$terms))
    
    p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, 
                                                        y = terms)) + ggplot2::geom_vline(xintercept = base::c(1:base::length(base::unique(ggplot_df$cluster))), 
                                                                                          color = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                   "cluster"])), alpha.f = 0.5), size = 1.5) + ggplot2::geom_point(ggplot2::aes(fill = cluster), color = "white", stroke = 1.5, size = 6, pch = 21) + ggplot2::scale_fill_manual(values = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                    "cluster"]))) + ggplot2::scale_color_manual(values = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "cluster"])))) + ggplot2::scale_x_discrete(limits = base::as.character(base::levels(ggplot_df$cluster))) + ggplot2::scale_y_discrete(limits = rev) + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(size = 0.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              color = "lightgrey")) + ggplot2::xlab("") + ggplot2::ylab("") + 
      ggplot2::ggtitle("Reactome enrichment")
    m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
    mm <- reshape2::melt(m)
    p2 <- ggplot2::ggplot(mm, ggplot2::aes(Var1, Var2, fill = value)) + 
      ggplot2::geom_tile(color = "white", height = 1) + ggplot2::scale_fill_gradientn(colors = (grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, 
                                                                                                                                                               name = "RdBu"))))(base::length(base::seq(-2, 2, by = 0.1)))) + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::theme_light() + 
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      ggplot2::scale_x_discrete(limits = base::as.character(base::levels(ggplot_df$cluster)))
  }else{      
    ggplot_df <- base::data.frame(terms = base::as.vector(base::as.matrix(top_GO)), 
                                  cluster = base::rep(base::colnames(top_GO), each = base::nrow(top_GO)), 
                                  val = base::factor(base::rep(base::colnames(top_GO), 
                                                               each = base::nrow(top_GO))) %>% base::as.numeric())
    ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), 
    ]
    
    
    
    p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, 
                                                        y = terms)) + ggplot2::geom_vline(xintercept = base::c(1:base::length(base::unique(ggplot_df$cluster))), 
                                                                                          color = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                   "cluster"])), alpha.f = 0.5), size = 1.5) + ggplot2::geom_point(ggplot2::aes(fill = cluster), color = "white", stroke = 1.5, size = 6, pch = 21) + ggplot2::scale_fill_manual(values = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                    "cluster"]))) + ggplot2::scale_color_manual(values = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "cluster"])))) + ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "cluster"]))) + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(size = 0.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       color = "lightgrey")) + ggplot2::xlab("") + ggplot2::ylab("") + 
      ggplot2::ggtitle("Reactome enrichment")
    m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
    mm <- reshape2::melt(m)
    p2 <- ggplot2::ggplot(mm, ggplot2::aes(Var1, Var2, fill = value)) + 
      ggplot2::geom_tile(color = "white", height = 1) + ggplot2::scale_fill_gradientn(colors = (grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, 
                                                                                                                                                               name = "RdBu"))))(base::length(base::seq(-2, 2, by = 0.1)))) + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::theme_light() + 
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), 
                                                                                   "cluster"])))
  }
  cp <- egg::ggarrange(p, p2, ncol = 1, heights = c(3, 1))
  Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], 
                                      hcobject[["global_settings"]][["save_folder"]], "/DJplot_Reactome_top_", 
                                      top, ".pdf"), width = 12, height = 12)
  print(cp)
  grDevices::dev.off()
  openxlsx::write.xlsx(x = res, file = base::paste0(hcobject[["working_directory"]][["dir_output"]], 
                                                    hcobject[["global_settings"]][["save_folder"]], "/reactome_enrichment.xlsx"), 
                       overwrite = T)
  output <- list()
  output[["p"]] <- cp
  output[["result"]] <- top_GO
  output[["enrichment"]] <- res
  hcobject[["integrated_output"]][["enrichments"]][["top_Reactome"]] <<- output
}


}