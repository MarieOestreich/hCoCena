#' Gene To Cluster Dictionary
#' 
#' The function maps the gene names to their corresponding cluster.
#' @return A data frame with two columns, the first containing gene names as strings, the second containing cluster colours as strings.
#' @export

GeneToCluster <- function(cluster_information = hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]){
  gtc <- base::do.call(rbind, base::apply(cluster_information,1, function(x){
    tmp <- x["gene_n"] %>%
      base::strsplit(., split = ",") %>%
      base::unlist(.)
    base::data.frame(gene = tmp, color = base::rep(x["color"], base::length(tmp)))
  }))
  return(gtc)
}