#' Import Clusters From File
#' 
#' Uses an imported clustering model instead of clustering the integrated network.
#' 	The model must be saved as two columns, the first containing genes, the second cluster colors (this is the format export_clusters() exports to).
#' @param file File path.
#' @param sep The separator of the file. Default is tab-separated.
#' @param header A Boolean. Whether or not the file has headers (column names).
#' @export

import_clusters <- function(file, sep = "\t", header = T){
  
  
  gtc <- readr::read_delim(file = file , delim = sep, col_names = header)
  base::colnames(gtc) <- base::c("gene", "cluster")
  gtc[] <- base::lapply(gtc, base::as.character)
  
  new_cluster_info <- NULL
  for(c in base::unique(gtc$cluster)){
    tmp <- gtc[gtc$cluster == c, ]
    gene_no <- base::nrow(tmp)
    gene_n <- base::paste0(tmp$gene, collapse = ",")
    if(c == "white"){
      cluster_included <- "no"
      vertexsize <- 1
    }else{
      cluster_included <- "yes"
      vertexsize <- 3
    }
    color <- c
    conditions <- base::paste0(base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])[1:(base::ncol(hcobject[["integrated_output"]][["GFC_all_layers"]])-1)], collapse = "#")
    gfc_means = hcobject[["integrated_output"]][["GFC_all_layers"]][hcobject[["integrated_output"]][["GFC_all_layers"]][["Gene"]] %in% tmp$gene,] %>%
      dplyr::select(-Gene) %>%
      base::colMeans()
    grp_means = base::paste0(base::round(gfc_means,3) , collapse = ",")
    
    new_cluster_info <- base::rbind(new_cluster_info,
                              base::data.frame(clusters = c,
                                         gene_no = gene_no,
                                         gene_n = gene_n,
                                         cluster_included = cluster_included,
                                         color = c,
                                         conditions = conditions,
                                         grp_means = grp_means,
                                         vertexsize = vertexsize,
                                         stringsAsFactors = F))
    
  }
  hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]] <<- new_cluster_info
  
}