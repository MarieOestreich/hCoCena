#' Visualize Gene Expression
#' 
#' Plots the mean expression values per condition for the given genes for each dataset as a heatmap. Values are scaled across rows.
#' @param genes A vector of strings giving the genes to be plotted.
#' @param name A string giving the name for the save file WITHOUT file ending.
#' @param width Width of the PDF, default ist 15, only change if plots overlap in PDF.
#' @param height Height of the PDF, default ist 10, only change if plots overlap in PDF.
#' @export

visualize_gene_expression <- function(genes, name = NULL, width = 15, height = 10){
  
  gtc <- GeneToCluster() %>% dplyr::filter(., gene %in% genes) %>% dplyr::filter(., !color == "white")
  
  if(nrow(gtc) == 0){
    message("None of these genes are present in the network.")
    stop()
  }else{
    message(base::nrow(gtc), " out of ", base::length(genes), " requested genes are present in the network.")
  }
  
  genes <- gtc$gene
  
  
  plotls <- list()
  
  for(x in base::seq_along(hcobject[["layers"]])){
    
    exp <- hcobject[["data"]][[base::paste0("set",x, "_counts")]][genes, ]
    
    
    hm_anno <- hcobject[["data"]][[base::paste0("set",x, "_anno")]][base::colnames(exp),] %>% 
      dplyr::select(., hcobject[["global_settings"]][["voi"]]) 
    base::colnames(hm_anno) <- "voi"
    
    #cluster_hm_colnames <- base::colnames(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]]@ht_list[[1]]@matrix)
    #cluster_hm_colnames <- cluster_hm_colnames[ComplexHeatmap::column_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]])]
    #cluster_hm_colnames <- cluster_hm_colnames[cluster_hm_colnames %in% hm_anno[["voi"]]]

    conditions <- base::unique(hm_anno$voi) %>% base::sort(.)
    
    mexp <- base::lapply(conditions, function(y){
      samples <- dplyr::filter(hm_anno, voi == y) %>% base::rownames()
      tmp_exp <- dplyr::select(exp, samples)
      tmp_mexp <- base::data.frame(V1 = base::apply(tmp_exp, 1, base::mean))
      base::colnames(tmp_mexp) <- y
      return(tmp_mexp)
    })%>% rlist::list.cbind()

    
    base::rownames(mexp) <- base::lapply(base::rownames(mexp), function(r){
      color <- dplyr::filter(gtc, gene == r) %>% dplyr::pull(., "color")
      return(paste0(r, " [", color, "]"))
    }) %>% base::unlist()
    
    hm <- pheatmap::pheatmap(mexp, 
                              scale = "row", 
                              cluster_rows = F,
                              cluster_cols = F, 
                              cellwidth = 60, 
                              cellheight = 30, 
                              main = base::paste0(hcobject[["layers_names"]][x], ": mean expression per condition"), 
                              color = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG")))(51))
    plotls[[x]] <- hm[["gtable"]]
  }
  
  cp <- cowplot::plot_grid(plotlist = plotls, ncol = base::length(hcobject[["layers"]]))
  cp

  if(base::is.null(name)){
    message("Cannot save the file since no unique file name was provided (See function parameter 'name').")
  }else{
    Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", name,
                                ".pdf"), width = width, height = height)

  #graphics::plot(cp)
  print(cp)

  grDevices::dev.off()
  }

}