#' Transfiction Factor Overrepresentation Analysis
#' 
#' Runs a ChEA3 transcription factor enrichment analysis for every gene module detected in the network. 
#' 	It filters the ranked enriched transcription factors for the 5 highest ranking ones and their 5 highest ranking targets.
#' 	MeanRank was chosen as a ranking method.
#' @export

TF_overrep <- function(){

  topTF = 5
  topTarget = 5
  output <- list()
  gtc <- GeneToCluster() 
  base::colnames(gtc) <- base::c("gene", "cluster")
  for(c in base::unique(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]][["color"]])){
    if(!c=="white"){
      genes = dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], color == c)%>%
        dplyr::pull(., "gene_n")%>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)
      
      url = "https://maayanlab.cloud/chea3/api/enrich/"  #"https://amp.pharm.mssm.edu/chea3/api/enrich/"
      encode = "json"
      payload = list(query_name = "myQuery", gene_set = genes)
      
      #POST to ChEA3 server
      response = httr::POST(url = url, body = payload, encode = encode)
      json = httr::content(response, as = "text")
      
      #results as list of R dataframes
      results = jsonlite::fromJSON(json)
      results <- results$`Integrated--meanRank`
      gtc_not_white <- gtc[!gtc$cluster == "white",]
      results <- dplyr::filter(results, TF %in% gtc_not_white$gene)
      # extract those from meanRank since meanRank scored as best method:
      resultlist <- list()
      for(i in 1:topTF){
        if(i > base::length(results$TF)){
          next
        }else{
          tf <- results$TF[i]
          overlapping_genes <- results$Overlapping_Genes[i]%>%
            base::strsplit(., split = ",")%>%
            base::unlist(.)
          resultlist[[tf]] <- list(TF = tf, targets = overlapping_genes[1:topTarget])
        }
      }
      
      # filter for only those that are present in our network:
      output[[c]] <- resultlist
    }
  }

  hcobject[["integrated_output"]][["TF_overrep_results"]] <<- output
}



#' Network-Wide Transcription Factor Overrepresentation Analysis
#' 
#' Returns the transcription factors that have the most enriched targets network-wide, including their targets.
#' @param topTF The number of transcription factors with the highest number of enricht targets in the network. Default is 100.
#' @param topTarget Per transcription factor the number of top most enriched targets to return. Default is 30.
#' @export


TF_enrich_all <- function(topTF = 100, 
							topTarget = 30){

  gtc <- GeneToCluster() 
  base::colnames(gtc) <- base::c("gene", "cluster")
  genes <- gtc$gene
  
  
  url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
  encode = "json"
  payload = list(query_name = "myQuery", gene_set = genes)
  
  #POST to ChEA3 server
  response = httr::POST(url = url, body = payload, encode = encode)
  json = httr::content(response, as = "text")
  
  #results as list of R dataframes
  results = jsonlite::fromJSON(json)
  results <- results$`Integrated--meanRank`
  gtc_not_white <- gtc[!gtc$cluster == "white",]
  results <- dplyr::filter(results, TF %in% gtc_not_white$gene)
  # extract those from meanRank since meanRank scored as best method:
  resultlist <- list()
  for(i in 1:topTF){
    
    tf <- results$TF[i]
    overlapping_genes <- results$Overlapping_Genes[i]%>%
      base::strsplit(., split = ",")%>%
      base::unlist(.)
    # genes to which this TF has an edge:
    edgelist <- dplyr::filter(hcobject[["integrated_output"]][["combined_edgelist"]], V1 == tf | V2 == tf)
    edgelist <- base::c(edgelist[,1] %>% base::as.character(), edgelist[,2] %>% base::as.character())
    edgelist <- edgelist[!edgelist == base::as.character(i)]
    
    overlapping_genes <- overlapping_genes[overlapping_genes %in% edgelist]

    if(base::length(overlapping_genes) > topTarget){
      resultlist[[tf]] <- list(TF = tf, targets = overlapping_genes[1:topTarget])
    }else{
      resultlist[[tf]] <- list(TF = tf, targets = overlapping_genes)
    }
    
  }
  hcobject[["integrated_output"]][["enrichall"]] <<- resultlist
}


#' Visualize Transcription Factor Enrichment Results
#' 
#' Results are visualized as circular plots, one for each cluster. Each cell in the plot is labelled with a gene name. 
#' 	The names of the transcription factors are highlighted in turquoise, those of target genes are written in black. 
#' 	The colour of the cell symbolises in which module the target or the transcription factor can be found in. 
#' 	If a gene is a target of one the enriched transcription factors, a link is connecting their cells. 
#' 	If there is an edge in the constructed co-expression network connecting those two genes, that link is non-transparent. 
#' 	If the TF/target-pair was found by ChEA3 but is not represented by an edge in the network, the link is depicted as slightly transparent.
#' @export

plot_TF_enrichment <- function(){

  gtc <- GeneToCluster() 
  base::colnames(gtc) <- base::c("gene", "cluster")
  dflist <- list()
  exp_plot_list <- list()
  TFs <- NULL
  for(c in base::names(hcobject[["integrated_output"]][["TF_overrep_results"]])){
    tmp <- NULL
    exp_plot_df <- base::names(hcobject[["integrated_output"]][["TF_overrep_results"]][[c]])
    for(t in base::names(hcobject[["integrated_output"]][["TF_overrep_results"]][[c]])){
      TFs <- base::c(TFs, t)
      clt <- dplyr::filter(gtc, gene == t)%>%
        dplyr::pull(., "cluster")
      tmp_df <- base::data.frame(TF = base::rep(t, base::length(hcobject[["integrated_output"]][["TF_overrep_results"]][[c]][[t]][["targets"]])), 
                           Target = hcobject[["integrated_output"]][["TF_overrep_results"]][[c]][[t]][["targets"]],
                           ClusterTF = base::rep(clt, base::length(hcobject[["integrated_output"]][["TF_overrep_results"]][[c]][[t]][["targets"]])))
      
      base::colnames(tmp_df) <- base::c("TF", "Target", "ClusterTF")
      tmp <- base::rbind(tmp, tmp_df)
      exp_plot_df <- base::c(exp_plot_df, hcobject[["integrated_output"]][["TF_overrep_results"]][[c]][[t]][["targets"]])
    }
    exp_plot_df <- base::as.data.frame(base::unique(exp_plot_df))
    base::colnames(exp_plot_df) <- c
    dflist[[c]] <- tmp
    exp_plot_list[[c]] <- exp_plot_df
  }
  
  
  TFs <- base::unique(base::as.character(TFs))
  edgelist <- hcobject[["integrated_output"]][["combined_edgelist"]]
  edgelist$merged <- base::paste0(base::as.character(edgelist$V1), base::as.character(edgelist$V2))
  edgelist$merged2 <- base::paste0(base::as.character(edgelist$V2), base::as.character(edgelist$V1))
  
  Cairo::CairoPDF(file = paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/circos_plots.pdf"), width = 12, 
                    height = 7, onefile = T)

  # loop in pairs of two to create douple plots: 
  for(n in 1:base::length(dflist)){
    # i1 <- n*2-1
    # i2 <- n*2
    # set layout to plot two plots each in horizontal arrangement
    graphics::layout(base::matrix(1:2, 1, 2)) 
    exp_p <- plot_TF(exp_plot_list[[n]])
    ComplexHeatmap::draw(exp_p)
    # for(j in c(i1, i2)){
      # catching 'out-of-bounds':
      # if(j > length(dflist)){
      #   break 
      # }
      # create link data frame:
      fromto <- dflist[[n]]
      # only consider no TF targets:
      fromto <- dplyr::filter(fromto, !Target %in% TF)
      fromto <- fromto[stats::complete.cases(fromto),]
      fromto <- base::unique(fromto)
      
      
      # dataframe that associates each gene with its cluster colour:
      NodeToColor <- base::rbind(base::data.frame(gene = fromto$TF, color = fromto$ClusterTF),
                           base::data.frame(gene = fromto$Target, color = base::rep(base::names(dflist)[n], base::nrow(fromto)))) %>%
        base::unique()
      
      # create plot factors:
      factors <- base::unique(base::as.character(NodeToColor$gene))
      circlize::circos.par(points.overflow.warning=FALSE)
      circlize::circos.initialize(factors, xlim = c(0, 1)) 
      circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = base::as.character(NodeToColor$color),
                             bg.border = NA ) 
      # add sector labels:
      g <- circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x,y){
        xlim = circlize::get.cell.meta.data("xlim")
        ylim = circlize::get.cell.meta.data("ylim")
        sector.name = circlize::get.cell.meta.data("sector.index")
        if(sector.name %in% TFs){
          circlize::circos.text(base::mean(xlim), base::mean(ylim)+2.5, sector.name, facing = "inside", niceFacing = T, cex = .9, 
                      col = "turquoise3", font = 2)
        }else{
          circlize::circos.text(base::mean(xlim), base::mean(ylim)+2.5, sector.name, facing = "inside", niceFacing = T, cex = .9)
        }
        
      })
      # add links
      for(i in 1:base::nrow(fromto)) {
        merged <- base::paste0(base::as.character(fromto[i,1]), base::as.character(fromto[i,2]))
        if(merged %in% edgelist$merged | merged %in% edgelist $merged2){
          g <- circlize::circos.link(sector.index1 =  base::as.character(fromto[i,1]), c(0.45, 0.55),
                           sector.index2 =  base::as.character(fromto[i,2]), c(.92), 
                           col = base::as.character(fromto[i,3]),
                           directional = 1,
                           arr.width = .1,
                           arr.length = .1)
        }else{
          g <- circlize::circos.link(sector.index1 =  base::as.character(fromto[i,1]), c(0.45, 0.55),
                           sector.index2 =  base::as.character(fromto[i,2]), c(.92), 
                           col = makeTransparent(base::as.character(fromto[i,3]), alpha = 30),
                           directional = 1,
                           arr.width = .1,
                           arr.length = .1)
        }
        
        
      }
      graphics::title(base::names(dflist)[n])
      circlize::circos.clear()
    # }
      
  }
  grDevices::dev.off()
}



