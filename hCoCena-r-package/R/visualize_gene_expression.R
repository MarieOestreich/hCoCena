#' Visualize Gene Expression
#' 
#' Plots the mean expression values per condition for the given genes for each dataset as a heatmap. Values are scaled across rows.
#' @param genes A vector of strings giving the genes to be plotted.
#' @param save Bolean to determine whether the plot should be saved as PDF. TRUE by default.
#' @param name A string giving the plot title and the name for the save file (WITHOUT file ending).
#' @param width Width of the PDF, default ist 15, only change if plots overlap in PDF.
#' @param height Height of the PDF, default ist 10, only change if plots overlap in PDF.
#' @export

visualize_gene_expression <- function(genes, name = NULL, width = 15, height = 10, save=T){
  # Filter for genes present in the network
  
  gtc <- GeneToCluster() %>% dplyr::filter(., gene %in% genes) %>% dplyr::filter(., !color == "white")
  
  if(nrow(gtc) == 0){
    message("None of these genes are present in the network.")
    stop()
  }else{
    message(base::nrow(gtc), " out of ", base::length(genes), " requested genes are present in the network.")
  }
  genes <- gtc$gene
  
  
  # Heatmap
  
  plotls <- NULL
  cp <- NULL 
  
  for(x in base::seq_along(hcobject[["layers"]])){
    
    exp <- hcobject[["data"]][[base::paste0("set",x, "_counts")]][genes, ]
    
    hm_anno <- hcobject[["data"]][[base::paste0("set",x, "_anno")]][base::colnames(exp),] %>% 
      dplyr::select(., dplyr::all_of(hcobject[["global_settings"]][["voi"]])) 
    base::colnames(hm_anno) <- "voi"

    conditions <- base::unique(hm_anno$voi) %>% base::sort(.)
    
    mexp <- base::lapply(conditions, function(y){
      samples <- base::subset(hm_anno, voi == y) %>% base::rownames()
      tmp_exp <- dplyr::select(exp, dplyr::all_of(samples))
      tmp_mexp <- base::data.frame(V1 = base::apply(tmp_exp, 1, base::mean))
      base::colnames(tmp_mexp) <- y
      return(tmp_mexp)
    })%>% rlist::list.cbind()
    
    base::rownames(mexp) <- base::lapply(base::rownames(mexp), function(r){
      color <- dplyr::filter(gtc, gene == r) %>% dplyr::pull(., "color")
      return(paste0(r, " [", color, "]"))
    }) %>% base::unlist()

    hm <- ComplexHeatmap::pheatmap(mat = as.matrix(mexp), 
                                   scale = "row", 
                                   cluster_rows = F, 
                                   cluster_cols = F, 
                                   cellwidth = 30, 
                                   cellheight = 15, 
                                   angle_col = "90", 
                                   legend = if(x==1){F}else{T},
                                   heatmap_legend_param = list(title = "scaled mean expr."), 
                                   main = hcobject[["layers_names"]][x], 
                                   fontsize = 10, 
                                   color = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG")))(51), 
                                   treeheight_col = 25, treeheight_row = 25)
    plotls <- plotls + hm
  }
  
  heatmap_title <- stringr::str_replace_all(string = name, pattern = "_", replacement = " ")
  
  if(save){
    if(base::is.null(name)){
      message("Cannot save the file since no unique file name was provided (See function parameter 'name').")
    }else{
      Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", name, ".pdf"),
                      width = width, height = height)
      ComplexHeatmap::plot.HeatmapList(plotls, column_title = heatmap_title, column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"))
      grDevices::dev.off()
    }
  }
  
  graphics::plot.new()
  ComplexHeatmap::plot.HeatmapList(plotls, column_title = heatmap_title, column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"))
   
}
