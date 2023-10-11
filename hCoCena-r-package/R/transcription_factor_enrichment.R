#' Transcription Factor Over-representation Analysis
#' 
#' Runs a ChEA3 transcription factor (TF) enrichment analysis for every selected gene module. 
#' 	It filters the ranked enriched TFs for the top highest ranking ones and their top highest ranking targets.
#' 	MeanRank was chosen as a ranking method.
#' Results are visualized as circular plots, one for each module. 
#' 	The names of the TFs are highlighted in turquoise, those of target genes are written in black. 
#' 	The colour of the cell symbolises in which module the respective gene can be found in. 
#' 	If a gene is a target of one of the enriched transcription factors, a link is connecting their cells. 
#' 	If there is an edge in the constructed co-expression network connecting those two genes, that link is non-transparent. 
#' 	If the TF/target-pair was found by ChEA3 but is not represented by an edge in the network, the link is depicted as slightly transparent.
#' @param topTF Integer. The number of top ranking TFs to return per cluster. Default is 5.
#' @param topTarget Integer. The number of top ranking targets to return per TF. Default is 5.
#' @param clusters Either "all" (default) or a vector of clusters as strings. Defines for which clusters to perform the analysis.
#' @export

TF_overrep <- function(clusters = "all", topTF = 5, topTarget = 5){

  output <- list()
  gtc <- GeneToCluster() 
  base::colnames(gtc) <- base::c("gene", "cluster")
  
  all_clusters <- base::unique(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]][["color"]])
  all_clusters <- all_clusters[all_clusters != "white"]
  
  if (clusters[1] == "all") {
    clusters <- all_clusters
  }
  
  tt_list <- list()
  exp_plot_list <- list()
  TFs <- NULL
  
  for(c in clusters){
    genes = dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], color == c) %>%
      dplyr::pull(., "gene_n")%>%
      base::strsplit(., split = ",") %>%
      base::unlist(.)
    
    url = "https://maayanlab.cloud/chea3/api/enrich/"
    encode = "json"
    payload = list(query_name = "myQuery", gene_set = genes)
    
    #POST to ChEA3 server
    response = httr::POST(url = url, body = payload, encode = encode)
    json = httr::content(response, as = "text")
    
    #results as list of R dataframes
    results <- jsonlite::fromJSON(json)
    results <- results$`Integrated--meanRank`
    gtc_not_white <- gtc[!gtc$cluster == "white",]
    results <- dplyr::filter(results, TF %in% gtc_not_white$gene)
    
    # extract those from meanRank since meanRank scored as best method:
    resultlist <- list()
    for(i in 1:topTF){
      if(i > base::length(results$TF)){
        next
      }else{
        if(results$Overlapping_Genes[i] == ""){
          next
        }else{
          tf <- results$TF[i]
          overlapping_genes <- results$Overlapping_Genes[i]%>%
            base::strsplit(., split = ",")%>%
            base::unlist(.)
          resultlist[[tf]] <- list(TF = tf, targets = overlapping_genes[1:topTarget])
        }
      }
    }
    output[[c]] <- resultlist
    
    # prepare TF-target dataframe for visualization
    if (base::length(resultlist) == 0){
      next
    }
    
    tt_df <- NULL
    exp_plot_df <- base::names(resultlist)
    
    for(x in base::names(resultlist)){
      TFs <- base::c(TFs, x)
      clt <- dplyr::filter(gtc, gene == x) %>% dplyr::pull(., "cluster")
      tmp_df <- base::data.frame(TF        = base::rep(x, base::length(resultlist[[x]][["targets"]])), 
                                 Target    = resultlist[[x]][["targets"]],
                                 ClusterTF = base::rep(clt, base::length(resultlist[[x]][["targets"]])))
      tmp_df <- tmp_df[stats::complete.cases(tmp_df), ]
      tt_df <- base::rbind(tt_df, tmp_df)
      exp_plot_df <- base::c(exp_plot_df, resultlist[[x]][["targets"]])
    }
    
    exp_plot_df <- base::as.data.frame(base::unique(exp_plot_df))
    base::colnames(exp_plot_df) <- c
    
    tt_list[[c]] <- tt_df
    exp_plot_list[[c]] <- exp_plot_df
  }
  
  # Save TF enrichment results
  hcobject[["integrated_output"]][["TF_overrep_results"]] <<- output
  
  
  # Generate plots
  
  if(length(tt_list) == 0){
    stop("No transcription factors found to be enriched for any of the modules.")
  }
  
  TFs <- base::unique(base::as.character(TFs))
  edgelist <- hcobject[["integrated_output"]][["combined_edgelist"]]
  edgelist$merged  <- base::paste0(base::as.character(edgelist$V1), base::as.character(edgelist$V2))
  edgelist$merged2 <- base::paste0(base::as.character(edgelist$V2), base::as.character(edgelist$V1))
  
  Cairo::CairoPDF(file = paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/TF_overrep.pdf"), 
                  width = 12, height = 7, onefile = T)
  
  # loop in pairs of two to create double plots: 
  for(n in 1:base::length(tt_list)){
    # i1 <- n*2-1
    # i2 <- n*2
    # set layout to plot two plots each in horizontal arrangement
    graphics::layout(base::matrix(1:2, 1, 2)) 
    exp_p <- plot_TF(exp_plot_list[[n]])
    ComplexHeatmap::draw(exp_p)
    
    # for(j in c(i1, i2)){
    # catching 'out-of-bounds':
    # if(j > length(tt_list)){
    #   break 
    # }
    # create link data frame:
    fromto <- tt_list[[n]]
    # only consider no TF targets:
    fromto <- dplyr::filter(fromto, !Target %in% TF)
    fromto <- fromto[stats::complete.cases(fromto),]
    fromto <- base::unique(fromto)
    
    
    # dataframe that associates each gene with its cluster colour:
    NodeToColor <- base::rbind(base::data.frame(gene = fromto$TF, color = fromto$ClusterTF),
                               base::data.frame(gene = fromto$Target, color = base::rep(base::names(tt_list)[n], base::nrow(fromto)))) %>%
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
      if(merged %in% edgelist$merged | merged %in% edgelist$merged2){
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
    graphics::title(base::names(tt_list)[n])
    circlize::circos.clear()
    # }
    
  }
  grDevices::dev.off()
}


#' Network-Wide Transcription Factor Over-representation Analysis
#' 
#' Returns the transcription factors that have the most enriched targets network-wide, including their targets.
#' @param topTF The number of transcription factors with the highest number of enriched targets in the network. Default is 100.
#' @param topTarget Per transcription factor the number of top most enriched targets to return. Default is 30.
#' @export


TF_enrich_all <- function(topTF = 100, topTarget = 30){

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
  results <- jsonlite::fromJSON(json)
  results <- results$`Integrated--meanRank`
  gtc_not_white <- gtc[!gtc$cluster == "white",]
  results <- dplyr::filter(results, TF %in% gtc_not_white$gene)
  
  # extract those from meanRank since meanRank scored as best method:
  resultlist <- list()
  for(i in 1:topTF){
    if(i > length(results$TF)){
      break
    }
    tf <- results$TF[i]
    overlapping_genes <- results$Overlapping_Genes[i]%>%
      base::strsplit(., split = ",")%>%
      base::unlist(.)
    
    # genes to which this TF has an edge:
    edgelist <- dplyr::filter(hcobject[["integrated_output"]][["combined_edgelist"]], V1 == tf | V2 == tf)
    edgelist <- base::c(edgelist[,1] %>% base::as.character(), edgelist[,2] %>% base::as.character())
    edgelist <- edgelist[!edgelist == base::as.character(i)]
    
    message('the transcription factor ', tf, ' has ', length(overlapping_genes), ' targets. It has a co-expression above the cutoff with ', length(overlapping_genes[overlapping_genes %in% edgelist]), ' of these targets. The others will be discarded.')
    overlapping_genes <- overlapping_genes[overlapping_genes %in% edgelist]

    if(base::length(overlapping_genes) > topTarget){
      resultlist[[tf]] <- list(TF = tf, targets = overlapping_genes[1:topTarget])
    }else{
      resultlist[[tf]] <- list(TF = tf, targets = overlapping_genes)
    }
    
  }
  hcobject[["integrated_output"]][["enrichall"]] <<- resultlist
}


