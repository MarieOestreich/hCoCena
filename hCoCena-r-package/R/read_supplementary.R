#' Read Supplementary Data
#' 
#' Function to read and collect the supplementary files set in set_supp_files().
#' @export

read_supplementary <- function () 
{
  if(!is.null(Tf))
    hcobject[["supplementary_data"]][["TF"]] <<- utils::read.delim(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                hcobject[["supplement"]][1]), header = TRUE, check.names = F)
  if(!is.null(Hall))
    hcobject[["supplementary_data"]][["hallmark"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                              hcobject[["supplement"]][2]))
  if(!is.null(Go))
    hcobject[["supplementary_data"]][["go"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                        hcobject[["supplement"]][3]))
  if(!is.null(Kegg))
    hcobject[["supplementary_data"]][["Kegg"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                          hcobject[["supplement"]][4]))
  if(!is.null(Reactome))
    hcobject[["supplementary_data"]][["Reactome"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                              hcobject[["supplement"]][5]))
}