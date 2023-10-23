#' Change Grouping Parameter
#' 
#' The variable by which the samples are grouped and based on which the GFCs are calculated can be changed and the cluster heatmap will be replotted. 
#'  This is particularly useful in cases where different variables are potential candidates for driving the genes’ expression changes in the data 
#'  and an explorative approach is required to decide on the most suitable one.
#'  Note that any previously generated column annotation will not be plotted, since the grouping will change.
#'  If you eventually decide on another grouping variable, please run the analysis again entirely with the changed "voi" from the very beginning. 
#' @param group_by A string giving the new grouping variable. This must be a name of a column present in all annotation files.
#' @export


change_grouping_parameter <- function(group_by){
  
  # check if grouping variables are present:
  if(base::length(group_by) == 1){
    for(i in 1:base::length(hcobject[["layers"]])){
      if(!group_by %in% base::colnames(hcobject[["data"]][[base::paste0("set", i, "_anno")]])){
        stop("Grouping variable not present as column name in all annotation tables.")
      }
    }
  }else{
    stop("Please provide only one grouping variable that is a column name present in ALL annotation files.")
  }
 
  
  # to store temporary GFCs per layer:
  sep_GFCs <- list()
  tmp_data <- list()
  # iterate over data sets:  
  for(i in 1:base::length(hcobject[["layers"]])){
    
    meta_data <- hcobject[["data"]][[base::paste0("set", i, "_anno")]]
    
    meta_data$regrouped <- base::paste0(meta_data[[group_by]],"_", hcobject[["layers_names"]][i])
    
    # mark controls:
    old_control <- hcobject[["global_settings"]][["control"]]
    if(!hcobject[["global_settings"]][["control"]] == "none"){
      # temporarily change control to none, since controls in new grouping variable is not known:
      hcobject[["global_settings"]][["control"]] <<- "none"
      print("Regrouped GFCs will be calculated without reference to controls.")
    }
    
    sep_GFCs[[i]] <- GFC_calculation(info_dataset = meta_data, grouping_v = "regrouped", x = i)
    tmp_data[[base::paste0("set", i, "_anno")]] <- meta_data
  }
  
  sep_GFCs <- purrr::reduce(sep_GFCs, dplyr::full_join, by = "Gene")
  sep_GFCs[base::is.na(sep_GFCs)] <- -hcobject[["global_settings"]][["range_GFC"]]
  gene_col <- sep_GFCs$Gene
  sep_GFCs$Gene <- NULL
  sep_GFCs$Gene <- gene_col
  
  
  # numerical annotation:
  
    # to be added
  
  # categorical annotation:
  
    # to be added

  hm <- replot_cluster_heatmap(GFCs = sep_GFCs, 
                               cluster_columns = T, 
                               cluster_rows = T, 
                               return_HM = T,
                               group = "regrouped",
                               file_name = paste0("Heatmap_modules_regrouped_", group_by, ".pdf"),
                               data = tmp_data)
  # reset to old control:
  hcobject[["global_settings"]][["control"]] <<- old_control
  

}