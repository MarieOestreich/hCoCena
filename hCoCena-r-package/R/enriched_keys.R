############ NOT GONE THROUGH YET ############
enriched_keys <- function(file = F, path, enrichment_keys = NULL,
                          cluster_info = integrated_output$cluster_calc$cluster_information, db = "GO"){
  output <- list()
    clusters <- unique(cluster_info$color)
    clusters <- clusters[!clusters == "white"]
    # if file == F, the enrichment_keys must be provided and are checked against GO enrichment terms
    if(!file){
      if(is.null(enrichment_keys)){
        print("file is set to FALSE but no enrichment keys were provided")
        return()
      }
      enriched_list <- list()
      categories_per_cluster <- NULL
      for(c in clusters){
        
        genes <- dplyr::filter(cluster_info, color == c)%>%
          dplyr::pull(., "gene_n")%>%
          base::strsplit(., split = ",")%>%
          BiocGenerics::unlist(.)
        
        if(db == "GO"){
          if(global_settings$organism %in% c("human", "Human")){
            enrich <- clusterProfiler::enrichGO(genes,
                                                OrgDb = "org.Hs.eg.db",
                                                keyType = "SYMBOL",
                                                ont = "BP",
                                                pvalueCutoff = 0.05)
          }
          if(global_settings$organism %in% c("mouse", "Mouse")){
            enrich <- clusterProfiler::enrichGO(genes,
                                                OrgDb = "org.Mm.eg.db",
                                                keyType = "SYMBOL",
                                                ont = "BP",
                                                pvalueCutoff = 0.05)
            
          }
          
        }else if(db == "KEGG"){
          if(global_settings$organism %in% c("human", "Human")){
            entrez <- clusterProfiler:: bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = F)$ENTREZID
            enrich <- clusterProfiler::enrichKEGG(entrez,
                                                  organism = "hsa",
                                                  pvalueCutoff = 0.05)
          }
          if(global_settings$organism %in% c("mouse", "Mouse")){
            entrez <- clusterProfiler:: bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = F)$ENTREZID
            enrich <- clusterProfiler::enrichKEGG(entrez,
                                                  organism = "mmu",
                                                  pvalueCutoff = 0.05)
          }
          
        }else{
          print("invalid database")
          return(NULL)
        }
        
        
        
        
        
        
        enriched_list[[c]] <- enrich

        
        cell_enrich <- list(counts = list(), genes = list())

        hits <- 0
        top_10 <- enrich@result[order(enrich@result$Count, decreasing = T),][1:20,]
        for(type in enrichment_keys){
          
          tmp <- dplyr::filter(top_10, grepl(type, top_10$Description,  ignore.case = T) == T) %>%
            dplyr::pull(., "geneID") %>%
            paste0(., collapse = "/") %>% 
            base::strsplit(., split = "/") %>% 
            unlist(.) %>% 
            unique(.)
          hits <- hits + length(tmp)
          cell_enrich$counts[[type]] <- length(tmp)
          cell_enrich$genes[[type]] <- tmp
          
        }
        tmp <- data.frame(matrix(unlist(cell_enrich$counts), ncol= length(cell_enrich$counts), byrow=T) %>% t(),stringsAsFactors=FALSE)%>%
          cbind(., names(cell_enrich$counts))%>%
          cbind(., rep(c, length(cell_enrich$counts)))
        
        sink(file = paste0(global_settings$save_folder, "/enrichedInKeys_",c,".txt"))
        for (n in names(cell_enrich$genes)){
          cat(NULL, sep = "\n")
          cat(n, sep = "\n\n")
          cat(cell_enrich$genes[[n]], sep = "\n")
        }
        
        sink()
        colnames(tmp) <- c("count", "cell_type", "cluster")

        
        tmp$hits <- rep(hits, nrow(tmp))
        
        
        categories_per_cluster <- rbind(categories_per_cluster, tmp)
      }
      output[["enrichlist"]] <- enriched_list
    }else{
      f <- read.csv(path)
      enrichment_keys <- colnames(f)
      categories_per_cluster <- NULL

      for(c in clusters){
        
        genes <- dplyr::filter(cluster_info, color == c)%>%
          dplyr::pull(., "gene_n")%>%
          base::strsplit(., split = ",")%>%
          BiocGenerics::unlist(.)
        
        cell_enrich <- list(counts = list(), genes = list())
        hits <- 0
        
        for (type in enrichment_keys){
          tmp <- genes[genes %in% f[,c(type)]]
          hits <- hits + length(tmp)
          cell_enrich$counts[[type]] <- length(tmp)
          cell_enrich$genes[[type]] <- tmp
        }
        
        tmp <- data.frame(matrix(unlist(cell_enrich$counts), ncol= length(cell_enrich$counts), byrow=T) %>% t(),stringsAsFactors=FALSE)%>%
          cbind(., names(cell_enrich$counts))%>%
          cbind(., rep(c, length(cell_enrich$counts)))
        
        sink(file = paste0(working_directory$dir_output, global_settings$save_folder, "/enrichedInKeys_",c,".txt"))
        for (n in names(cell_enrich$genes)){
          cat(NULL, sep = "\n")
          cat(n, sep = "\n\n")
          cat(cell_enrich$genes[[n]], sep = "\n")
        }
        
        sink()
        colnames(tmp) <- c("count", "cell_type", "cluster")
        if(hits == 0){
          tmp$count <- 0
        }else{
          tmp$count <- (tmp$count/hits)*100
        }
        
        tmp$hits <- rep(hits, nrow(tmp))
        
        
        categories_per_cluster <- rbind(categories_per_cluster, tmp)
        
      }
    }
    output[["categories_per_cluster"]] <- categories_per_cluster

    if("enriched_per_cluster" %in% names(hcobject[["satellite_outputs"]])){
    	hcobject[["satellite_outputs"]][["enriched_per_cluster2"]] <<- output
    }else{
    	hcobject[["satellite_outputs"]][["enriched_per_cluster"]] <<- output
    }
}