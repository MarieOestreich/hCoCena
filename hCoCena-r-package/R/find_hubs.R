#' Find Hub Genes
#' 
#' Hub genes are determined per cluster using a combined ranking based on weighted degree centrality, weighted closeness centrality and weighted betweenness centrality.
#'  A table of hub genes per cluster is returned as an excel file and a heatmap of the hub genes expression values is plotted per cluster and per dataset.
#' @param top An integer. All genes are ranked based on their hub potential, this parameter defines the number of the top ranked genes to be considered a hub gene. Default value is 10. 
#' @param save A Boolean. Whether or not a labelled hub-network per cluster is to be save to PDF. Default is FASLE.
#' @param tree_layout A Boolean. Whether or not to depict the network witht tree layout, implying a sort of hierarchical structure to the network. 
#' @param TF_only Either FALSE (default, all genes in clsuter are considered for hub genes), or "all" (all genes from transcriptionfactor supplementary file are considered for hub genes),
#'  or any gene category listed in the last column of the provided transcriptionfactor supplementary file (only that subgroup condired for hub genes).
#' @param Plot A Boolean. Wheather or not to plot the network (per cluster) with highlighted hub nodes. Default is FALSE.
#' @param clusters Either "all" (default) or a vector of cluster colours for which the hub detection should be performed.
#' @export

find_hubs <- function(clusters = c("all"),
                      top = 10, 
                      save = F, 
                      tree_layout = F, 
                      TF_only = F, 
                      Plot = F){
  
  gtc <- GeneToCluster()
  if(clusters[1] == "all"){
    clusters <- base::unique(gtc$color[!gtc$color == "white"])
  }

  hubs <- base::lapply(clusters, function(x){

    tmp <- hub_node_detection(cluster = x, top = top, save = save, tree_layout = tree_layout, TF_only = TF_only, Plot = Plot)
    return(tmp$hub_nodes)

  })
  
  
  hubs_df <- base::lapply(hubs, function(x){
    if(base::is.null(x)){
      x <- base::rep(" ", top)
    }
    if( base::length(x) < top ){
      x <- base::c(x, base::rep(" ", top - base::length(x)))
    }
    return(x)
  }) %>% rlist::list.cbind() %>% base::as.data.frame()

  base::colnames(hubs_df) <- clusters
  
  hubs_df[base::is.na(hubs_df)] <- " "

  for(col in base::colnames(hubs_df)){
    if(base::file.exists(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/hubgenes.xlsx"))){
      wb <- openxlsx::loadWorkbook(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/hubgenes.xlsx"))
      if(col %in% wb$sheet_names) { openxlsx::removeWorksheet(wb, col) }
      openxlsx::addWorksheet(wb, col)
      openxlsx::writeData(wb, sheet = col, dplyr::select(hubs_df, tidyselect::all_of(col)), colNames = T)
      openxlsx::saveWorkbook(wb,base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/hubgenes.xlsx"),overwrite = T)
    }else{
      openxlsx::write.xlsx(dplyr::select(hubs_df, tidyselect::all_of(col)), 
                       file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/hubgenes.xlsx"), 
                       sheetName = col, 
                       colNames = TRUE, rowNames = F, append = FALSE, overwrite = T)
    }

    if("hub_out" %in% base::names(hcobject[["satellite_outputs"]])){
      hcobject[["satellite_outputs"]][["hub_out"]][[col]] <<- dplyr::select(hubs_df, tidyselect::all_of(col))
    }else{
      hcobject[["satellite_outputs"]][["hub_out"]] <<- list()
      hcobject[["satellite_outputs"]][["hub_out"]][[col]] <<- dplyr::select(hubs_df, tidyselect::all_of(col))
    }
  }
  
  #plot_hub_exp(hubs_df = hubs_df)
  for(col in base::colnames(hubs_df)){
    if(!all(dplyr::pull(hubs_df, col)== " ")){
      visualize_gene_expression(genes = dplyr::pull(hubs_df, col), name = base::paste0(col, "_hubs"), width = 25)
    }
  }
}


