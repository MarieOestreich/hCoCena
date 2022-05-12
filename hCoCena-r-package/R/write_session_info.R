#' Write Session Information
#' 
#' The parameters of the analysis session are written to a text file to enhance reproducibility without keeping adn sharing a markdown for every analysis. 
#' It documents the name of the files and their location used as count and annotation files, the global settings set in the session, 
#' 	the layer settings set for each dataset as well as the cut-offs and the clustering algorithm used.
#' @export

write_session_info <- function(){

  f <- base::file(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/session_info.txt"), "w")
  # number of layers:
  base::writeLines(base::paste0("Number of layers: ", base::length(hcobject[["layers"]])), f, sep = "\n")
  # file names of data sets
  base::writeLines(base::c("Count files:"), f, sep = "\n")
  for(i in 1:base::length(hcobject[["layers"]])){
    if(hcobject[["working_directory"]][["dir_count_data"]] == FALSE){
      filepath <- base::paste0("    ", "-count data loaded from R object, not from file-")
    }else{
      filepath <- base::paste0("    ", hcobject[["working_directory"]][["dir_count_data"]], hcobject[["layers"]][[i]][[1]])
    }
    base::writeLines(base::c(filepath), f, sep = "\n")
  }
  base::writeLines(base::c("Annotation files:"), f, sep = "\n")
  for(i in 1:base::length(hcobject[["layers"]])){
    if(hcobject[["working_directory"]][["dir_annotation"]] == FALSE){
      filepath <- base::paste0("    ", "-annotation data loaded from R object, not from file-")
    }else{
      filepath <- base::paste0("    ", hcobject[["working_directory"]][["dir_annotation"]], hcobject[["layers"]][[i]][[2]])
    }
    
    base::writeLines(base::c(filepath), f, sep = "\n")
  }
  
  # global settings:
  base::writeLines(base::c("Global settings:"), f, sep = "\n")
  for(i in 1:base::length(hcobject[["global_settings"]])){
    if(!base::names(hcobject[["global_settings"]])[i] == "chosen_clustering_algo"){
      base::writeLines(base::paste0("    ", base::names(hcobject[["global_settings"]])[i], ": ", hcobject[["global_settings"]][[i]]), f, sep = "\n")
    }
  }
  
  # layer settings:
  base::writeLines(base::c("Layer settings:"), f, sep = "\n")
  for(i in 1:base::length(hcobject[["layers"]])){
    base::writeLines(base::paste0("    layer: ", i), f, sep = "\n")
    for(j in 1:base::length(hcobject[["layer_settings"]][[i]])){
      base::writeLines(base::paste0("        ", base::names(hcobject[["layer_settings"]][[i]])[j], ": ", hcobject[["layer_settings"]][[i]][[j]]), f, sep = "\n")
    }
  }
  
  # cutoffs:
  base::writeLines(base::c("Cut-offs used:"), f, sep = "\n")
  for(i in 1:base::length(hcobject[["layers"]])){
    base::writeLines(base::paste0("    layer ", i, ": ", hcobject[["cutoff_vec"]][i]), f, sep = "\n")
  }
  
  # clustering algo
  base::writeLines(base::c("Clustering algorithm used: "), f, sep = "\n")
  base::writeLines(base::paste0("    ", hcobject[["global_settings"]][["chosen_clustering_algo"]]), f, sep = "\n")
  base::close(f)
}