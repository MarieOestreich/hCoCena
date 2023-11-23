#' User Specific Cluster Profiling
#' 
#' Offers a user-defined enrichment analysis for the identified modules. There are two options:
#' 	1) The user provides their own enrichment file. 
#' 		This file should be provided in .csv format where the column names represent the categories (a.k.a. "keys", e.g. cell types), and the columns contain genes representative of that key. 
#' 		The genes do not have to be exclusive with regard to the key they belong to and the columns are not required to be of the same length. 
#' 		If such a file is provided, the "from_file" parameter needs to be set to TRUE and the path parameter needs to be set to the path at which to find the file. 
#' 	2) Instead of providing an enrichment file, the user can also choose a database - so far Gene Ontology and KEGG are supported – and define a vector of keys as strings. 
#' 		The function then conducts an enrichment analysis using clusterProfiler functions (Guangchuang Yu, 2012) and filters the enrichment results for terms including the defined keys. 
#' 		In this case, the "from_file" parameter needs to be set to FALSE and the "enrichment_keys" parameter must be set to the vector of keys, the modules should be screened for. 
#' 		By setting the parameter 'db = "GO"' or 'db = "KEGG"', the GO or KEGG database is chosen accordingly. 
#' 		Up to 2 user-defined enrichments are possible so far and they are saved as "enriched_per_cluster" and "enriched_per_cluster2" in hcobject$satellite_outputs. 
#' 		The enrichment results will be visualized as stacked bar plot annotations on the module heatmap when rerunning plot_cluster_heatmap().
#' @param from_file A Boolean. If the enrichment is based on option 1 (see above), set to TRUE.
#' @param path The path to the enrichment file. Can be ignored if 'from_file' is FALSE.
#' @param enrichment_keys A vector of keys, the database terms should be scanned for for each cluster. Will be ignored if 'from_file' is TRUE.
#' @param db "Go" to use Gene Ontology database, "Kegg" to use KEGG database, "Hallmark" to use Hallmark genesets or custom database. Will be ignored when 'from_file' is TRUE.
#' @param padj Method to use for multiple testing correction. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".  Default is "BH" (Benjamini-Hochberg).
#' @param qval q-value cutoff to define terms to consider after pathway enrichment.
#' @export 

user_specific_cluster_profiling <- function(from_file = F, 
                              							path = NULL, 
                              							enrichment_keys = NULL,
                              							db = c("Go", "Kegg", "Hallmark"),
                              							padj = "BH",
                              							qval = 0.1){

  	cluster_info <- hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]
  	output <- list()
    clusters <- base::unique(cluster_info$color)
    clusters <- clusters[!clusters == "white"]

    # if file == F, the enrichment_keys must be provided and are checked against GO enrichment terms
    if(!from_file){
      if(base::is.null(enrichment_keys)){
        print("file is set to FALSE but no enrichment keys were provided")
        return()
      }
      
      enriched_list <- list()
      categories_per_cluster <- NULL
      
      for(c in clusters){
        
        genes <- dplyr::filter(cluster_info, color == c) %>% dplyr::pull(., "gene_n") %>% base::strsplit(., split = ",") %>% base::unlist(.)
        
        if(stringr::str_to_title(db) %in% names(hcobject[["supplementary_data"]])){
            enrich <- clusterProfiler::enricher(genes,
                                                TERM2GENE = hcobject[["supplementary_data"]][[stringr::str_to_title(db)]],
                                                pAdjustMethod = padj,
                                                pvalueCutoff = 0.05,
                                                qvalueCutoff = qval)
            
        }else{
          print("invalid database")
          return(NULL)
        }
        
        enriched_list[[c]] <- enrich
        
        cell_enrich <- list(counts = list(), genes = list())

        hits <- 0
        top_terms <- enrich@result[enrich@result$qvalue <= qval,] # [1:20,]
        
        for(type in enrichment_keys){
          tmp <- dplyr::filter(top_terms, base::grepl(type, top_terms$Description,  ignore.case = T) == T) %>%
            dplyr::pull(., "geneID") %>%
            base::paste0(., collapse = "/") %>% 
            base::strsplit(., split = "/") %>% 
            base::unlist(.) %>% 
            base::unique(.)
          hits <- hits + base::length(tmp)
          cell_enrich$counts[[type]] <- base::length(tmp)
          cell_enrich$genes[[type]] <- tmp
        }

        tmp <- base::data.frame(base::matrix(base::unlist(cell_enrich$counts), 
                                             ncol = base::length(cell_enrich$counts), 
                                             byrow = T) %>% 
                                  base::t(),stringsAsFactors = FALSE) %>%
          base::cbind(., base::names(cell_enrich$counts))%>%
          base::cbind(., base::rep(c, base::length(cell_enrich$counts)))
        
        base::colnames(tmp) <- base::c("count", "cell_type", "cluster")
    
        tmp$hits <- base::rep(hits, base::nrow(tmp))
        
        if(hits == 0){
          tmp$count <- 0
        }else{
          tmp$count <- (tmp$count/tmp$hits)*100
        }
        
        categories_per_cluster <- base::rbind(categories_per_cluster, tmp)
      }
      output[["enrichlist"]] <- enriched_list
      
    }else{
      
      if(base::is.null(path)){
        stop("File path must be specified using the 'path' variable.")
      }
      f <- utils::read.csv(path)
      enrichment_keys <- base::colnames(f)
      categories_per_cluster <- NULL

      for(c in clusters){
        
        genes <- dplyr::filter(cluster_info, color == c)%>%
          dplyr::pull(., "gene_n")%>%
          base::strsplit(., split = ",")%>%
          base::unlist(.)
        
        cell_enrich <- list(counts = list(), genes = list())
        hits <- 0
        
        for (type in enrichment_keys){
          tmp <- genes[genes %in% f[, base::c(type)]]
          hits <- hits + base::length(tmp)
          cell_enrich$counts[[type]] <- base::length(tmp)
          cell_enrich$genes[[type]] <- tmp
        }
        
        tmp <- base::data.frame(base::matrix(base::unlist(cell_enrich$counts), ncol = base::length(cell_enrich$counts), byrow=T) %>% base::t(),stringsAsFactors=FALSE)%>%
          base::cbind(., base::names(cell_enrich$counts))%>%
          base::cbind(., base::rep(c, base::length(cell_enrich$counts)))

        base::colnames(tmp) <- c("count", "cell_type", "cluster")
        if(hits == 0){
          tmp$count <- 0
        }else{
          tmp$count <- (tmp$count/hits)*100
        }
        
        tmp$hits <- base::rep(hits, base::nrow(tmp))
        
        
        categories_per_cluster <- base::rbind(categories_per_cluster, tmp)
        
      }
    }
    output[["categories_per_cluster"]] <- categories_per_cluster
    base::attr(output[["categories_per_cluster"]],"hidden") <- list(from_file, path, enrichment_keys, db, padj, qval)
    if("enriched_per_cluster" %in% base::names(hcobject[["satellite_outputs"]])){
      if(identical(base::attr(output[["categories_per_cluster"]], "hidden"), base::attr(hcobject[["satellite_outputs"]][["enriched_per_cluster"]][["categories_per_cluster"]], "hidden")))
        hcobject[["satellite_outputs"]][["enriched_per_cluster"]] <<- output
      else
    	  hcobject[["satellite_outputs"]][["enriched_per_cluster2"]] <<- output
    }else{
    	hcobject[["satellite_outputs"]][["enriched_per_cluster"]] <<- output
    }
}
