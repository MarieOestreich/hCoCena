#' Read Supplementary Data
#' 
#' Function to read and collect the supplementary files set in set_supp_files().
#' @export

read_supplementary <- function (){
  sapply(names(hcobject[["supplement"]]), function(x))  
    if(x == "Tf" & !is.null(hcobject[["supplement"]][["Tf"]])) 
      hcobject[["supplementary_data"]][["TF"]] <<- utils::read.delim(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                  hcobject[["supplement"]][["Tf"]]), header = TRUE, check.names = F)
    
    if(x == "Hallmark" & !is.null(hcobject[["supplement"]][["Hallmark"]]))
      hcobject[["supplementary_data"]][["Hallmark"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                                hcobject[["supplement"]][["Hallmark"]]))
    
    if(x == "Go" & !is.null(hcobject[["supplement"]][["Go"]]))
      hcobject[["supplementary_data"]][["Go"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                          hcobject[["supplement"]][["Go"]]))
    if(x == "Kegg" & !is.null(hcobject[["supplement"]][["Kegg"]]))
      hcobject[["supplementary_data"]][["Kegg"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                            hcobject[["supplement"]][["Kegg"]]))
    if(x == "Reactome" & !is.null(hcobject[["supplement"]][["Reactome"]]))
      hcobject[["supplementary_data"]][["Reactome"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                                hcobject[["supplement"]][["Reactome"]]))
    if(!(x %in% c("Tf", "Hallmark", "Go", "Kegg", "Reactome")))
      if(grepl(".csv", hcobject[["supplement"]][[x]]))
        hcobject[["supplementary_data"]][[stringr::str_to_title(x)]] <<- utils::read.csv(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                                      hcobject[["supplement"]][[x]]), header = T, stringAsFactors = F, quote ="")
      else if(grepl(".gmt", hcobject[["supplement"]][[x]]))
        hcobject[["supplementary_data"]][["Reactome"]] <<- clusterProfiler::read.gmt(base::paste0(hcobject[["working_directory"]][["dir_reference_files"]], 
                                                                                                  hcobject[["supplement"]][[x]]))
      else print(paste0("invalid input format of database: ", x, ". Valid inputs are: .csv and .gmt files!"))
}
