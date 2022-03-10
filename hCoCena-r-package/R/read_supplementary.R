#' Read Supplementary Data
#' 
#' Function to read and collect the supplementary files set in set_supp_files().
#' @export

read_supplementary <- function(){
  
  
  hcobject[["supplementary_data"]][["TF"]] <<- utils::read.delim(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], hcobject[["supplement"]][1]),
                             header=TRUE,
                             check.names=F)
  
  hcobject[["supplementary_data"]][["hallmark"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], hcobject[["supplement"]][2]))
  
  hcobject[["supplementary_data"]][["go"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], hcobject[["supplement"]][3]))
  
}