#' Read Supplementary Data
#' 
#' Function to read and collect the supplementary files set in set_supp_files().
#' @export

read_supplementary <- function (){
  
  if(!is.null(hcobject[["supplement"]][["Tf"]])) 
    hcobject[["supplementary_data"]][["TF"]] <<- utils::read.delim(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                hcobject[["supplement"]][["Tf"]]), header = TRUE, check.names = F)
  
  if(!is.null(hcobject[["supplement"]][["Hallmark"]]))
    hcobject[["supplementary_data"]][["Hallmark"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                              hcobject[["supplement"]][["Hallmark"]]))
  
  if(!is.null(hcobject[["supplement"]][["Go"]]))
    hcobject[["supplementary_data"]][["Go"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                        hcobject[["supplement"]][["Go"]]))
  if(!is.null(hcobject[["supplement"]][["Kegg"]]))
    hcobject[["supplementary_data"]][["Kegg"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                          hcobject[["supplement"]][["Kegg"]]))
  if(!is.null(hcobject[["supplement"]][["Reactome"]]))
    hcobject[["supplementary_data"]][["Reactome"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                              hcobject[["supplement"]][["Reactome"]]))
}
