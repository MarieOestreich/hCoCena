#' Export Clustering
#' 
#' Exports your clustering to a text file that can be imported into other hCoCena runs. 
#'  Format: Two columns, the first one named gene (containing the gene names) and the second one named color (containing the corresponding cluster color to which the gene belonged).
#' @export

export_clusters <- function(){
  gtc <- GeneToCluster()
  readr::write_delim(gtc, file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/gtc.txt"), 
                     delim = "\t", 
                     col_names = T)
}