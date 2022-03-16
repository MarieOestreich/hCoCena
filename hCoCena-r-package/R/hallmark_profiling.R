#' Hallmark Profiling
#' 
#' Returns a cluster-wise Hallmark Database enrichment.
#' @param top Integer. The number of most strongly enriched terms to return per cluster.
#' @param padj Method to use for multiple testing correction. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".  Default is "BH" (Benjamini-Hochberg).
#' @param clusters Either "all" (default) or a vector of clusters as strings. Defines for which clusters to perform the enrichment.
#' @param qval Upper threshold for the adjusted p-value. Default is 0.05.
#' @export

Hallmark_profiling <- function(top, padj = "BH", clusters = c("all"), qval = 0.05){
  
  cluster_info <- hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]

  # gene universe: union of all genes present in the count files:
  universe <- base::lapply(1:base::length(hcobject[["layers"]]), function(x){
    return(base::rownames(hcobject[["data"]][[base::paste0("set", x, "_counts")]]))
  }) %>% base::unlist() %>% base::unique()
  
  if(clusters[1]=="all"){
    clusters <- base::unique(cluster_info$color)
    clusters <- clusters[!clusters == "white"]
  }
  
  
  top_GO <- list()
  res <- list()
  for(c in clusters){
    
    genes <- dplyr::filter(cluster_info, color == c)%>%
      dplyr::pull(., "gene_n")%>%
      base::strsplit(., split = ",")%>%
      base::unlist(.)

    enrich <- clusterProfiler::enricher(genes,
                                        TERM2GENE=hcobject[["supplementary_data"]][["hallmark"]],
                                        qvalueCutoff = qval,
                                        pAdjustMethod = padj,
                                        universe = universe)
    
    if(!base::is.null(enrich)){
      tmp <- enrich@result
      tmp <- dplyr::filter(tmp, qvalue <= qval)
      if(base::nrow(tmp) == 0){
        next
      }else{
        tmp <- tmp[base::order(tmp$qvalue, decreasing = F),]
        if(base::nrow(tmp) < top){
          diff <- top-base::nrow(tmp)
          top_GO[[c]] <- c(tmp$Description, rep(NA, diff))
        }else{
          top_GO[[c]] <- tmp$Description[1:top]
        }
        top_GO[[c]] <- tmp$Description[1:top]
        res[[c]] <- enrich@result
      }
    }
    
  }
  
  top_GO <- dplyr::bind_rows(top_GO)
  if(base::nrow(top_GO) == 1){
    top_GO <- base::t(top_GO)
  }
  ggplot_df <- base::data.frame(terms = base::as.vector(base::as.matrix(top_GO)), cluster = base::rep(base::colnames(top_GO), each = base::nrow(top_GO)),
                     val = base::factor(base::rep(base::colnames(top_GO), each = base::nrow(top_GO))) %>% base::as.numeric())
  ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df),]
  p <- ggplot2::ggplot(data = ggplot_df, ggplot2::aes(x = val, y = terms))+
    ggplot2::geom_vline(xintercept = base::c(1:base::length(base::unique(ggplot_df$cluster))), color = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), "cluster"])), alpha.f = 0.5), size = 1.5)+
    ggplot2::geom_point(ggplot2::aes(fill = cluster, color = cluster),size = 6, pch = 22)+
    ggplot2::scale_fill_manual(values = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), "cluster"])))+
    ggplot2::scale_color_manual(values = grDevices::adjustcolor(base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), "cluster"]))))+
    ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), "cluster"])))+
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90), panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line( size=.1, color="lightgrey"))+
    ggplot2::xlab("")+
    ggplot2::ylab("")+
    ggplot2::ggtitle("HALLMARK enrichment")

  m <- hcobject$integrated_output$cluster_calc$heatmap_cluster@ht_list[[1]]@matrix
  mm <- reshape2::melt(m)
  
  p2 <- ggplot2::ggplot(mm, ggplot2::aes(Var1, Var2, fill = value))+
    ggplot2::geom_tile(color = "white", height = 1)+
    #ggplot2::coord_equal()+
    ggplot2::scale_fill_gradientn(colors = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(base::length(base::seq(-2, 2, by = .1))))+
    ggplot2::ylab("")+
    ggplot2::xlab("")+
    ggplot2::theme_light()+
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle=90))+
    ggplot2::scale_x_discrete(limits = base::as.character(base::unique(ggplot_df[base::order(ggplot_df$val), "cluster"])))

  cp <- egg::ggarrange(p, p2, ncol=1, heights = c(3,1))

  
  Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/DJplot_HALLMARK_top_",
                                top, ".pdf"),
                  width = 12, height = 12)

  #graphics::plot(cp)
  print(cp)

  grDevices::dev.off()

  openxlsx::write.xlsx(x = res, file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/hallmark_enrichment.xlsx"), overwrite = T)

  output <- list()
  output[["p"]] <- cp
  output[["result"]] <- top_GO
  output[["enrichment"]] <- res
  hcobject[["integrated_output"]][["enrichments"]][["top_HALL"]] <<- output
}