#' Function To Provide Cluster Colours
#' @noRd

get_cluster_colours <- function(){

  col_vec <- c("coral", "gold", "steelblue", "lightgreen", "turquoise", "plum", "maroon", "seagreen", "wheat", "slategray", "lightblue", 
                     "orchid", "darkgreen",  "darkorange", "darkgrey", "indianred", "pink", "sandybrown",   "khaki",  "darkblue", "cadetblue",
                     "greenyellow","cyan", "thistle", "darkmagenta",  "red", "blue", "green", "yellow", "brown", "black", "darkgoldenrod", 
                     "cornsilk", "firebrick", "deeppink", "dodgerblue", "lightpink", "midnightblue", "slategray", "aquamarine", "chocolate", 
                     "darkred", "navy", "olivedrab", "peachpuff", "tomato", "snow")
  return(col_vec)

}


#' Leiden Clustering
#' 
#' Applies the Leiden community detection to the network.
#' @param g Then network, an igraph object.
#' @param num_it The number of iteration the algorithm is supposed to run.
#' @noRd

leiden_clustering <- function(g, num_it){
  
  color.cluster <- get_cluster_colours()
  
  set.seed(168575)
  
  # run Leiden algorithm ion network:
  partition <- leidenAlg::leiden.community(graph = g, n.iterations = num_it)

  # extract found clusters and the genes belonging to them:
  clusters_df <- base::data.frame(cluster = base::as.numeric(partition$membership), gene = partition$names)

  # the algo start enumeration at 0, increase to 1 for easier comaptibility with R fucntions:
  clusters_df$cluster <- clusters_df$cluster + 1

  # get gene counts per cluster:
  cluster_frequencies <- base::table(clusters_df$cluster) %>% base::as.data.frame() 

  # detect clusters large enough to be kept:
  clusters_to_keep <- dplyr::filter(cluster_frequencies, Freq >= hcobject[["global_settings"]][["min_nodes_number_for_cluster"]]) %>% 
    dplyr::pull(., "Var1")
  
  # define white clusters (those that are too small to be kept):
  clusters_df_white <- dplyr::filter(clusters_df, !cluster %in% clusters_to_keep)

  # remove white clusters:
  clusters_df <- dplyr::filter(clusters_df, cluster %in% clusters_to_keep)

  # inform how many clusters and accordingly how many genes were lost due to insufficient cluster size:
  print(base::paste0(base::length(base::unique(clusters_df_white$cluster)), 
               " cluster/s was/were smaller than the set minimum cluster size and therefore discarded.",
               " This removes ", base::nrow(clusters_df_white), " genes."))


  out <- base::lapply(base::unique(clusters_df$cluster), function(x){

    # extract genes present in current cluster:
    genes <- dplyr::filter(clusters_df, cluster == x) %>% 
      dplyr::pull(., "gene")

    # cluster name:
    col_clusters <- base::paste0("cluster ", x)
    # number of genes in cluster:
    col_gene_no <- dplyr::filter(clusters_df, cluster == x) %>% 
      base::nrow()
    # comma separated list of all genes in the cluster:
    col_gene_n <- genes %>% base::paste0(., collapse = ",")
    # is the cluster included in the network (i.e., is it large enough):
    col_cluster_included <- "yes"
    # cluster color:
    col_color <- color.cluster[x]
    # order of conditions for following GFC values:
    col_conditions <- base::colnames(dplyr::select(hcobject[["integrated_output"]][["GFC_all_layers"]], -Gene)) %>% 
      base::paste0(., collapse = "#")
    # mean GFCs of the cluster genes per sample group:
    col_grp_means <- dplyr::filter(hcobject[["integrated_output"]][["GFC_all_layers"]], Gene %in% genes) %>% 
      dplyr::select(-Gene) %>% 
      base::apply(., 2, base::mean) %>% 
      base::round(., digits = 3) %>% 
      base::paste0(., collapse = ",")
    # collect all information:
    out <- base::data.frame(clusters = col_clusters,
                      gene_no = col_gene_no,
                      gene_n = col_gene_n,
                      cluster_included = col_cluster_included,
                      color = col_color,
                      conditions = col_conditions,
                      grp_means = col_grp_means,
                      vertexsize = 3)
    return(out)
    }) %>% rlist::list.rbind()
  return(out)
}

    
#' Internal Function Used In cluster_calculation()
#' @noRd

cluster_calculation_internal <- function(graph_obj, 
                                          algo, 
                                          case, 
                                          it = 1) {

    library(igraph)

    if(algo == "cluster_leiden"){
      cfg <- leiden_clustering(g = graph_obj, num_it = no_of_iterations)
    }
    cfg <- base::get(algo)(graph_obj)

    mod_score <- igraph::modularity(graph_obj, base::as.numeric(cfg$membership)) 

    mod_df <- base::data.frame(modularity_score = mod_score, 
                                cluster_algorithm = algo, 
                                stringsAsFactors = F)

    # making switch so that in the end when only the best algorithm is to be used then the same function can be used
    output <- base::switch(case, best = cfg$membership, test = mod_df, final = cfg)

    print(base::paste0(algo," algorithm tested"))
    return(output)
  }


#'Internal Function Used In cluster_calculation()
#' @noRd

gfc_mean_clustergene <- function(rownum, cluster_df, gfc_dat){

  d1 <- cluster_df[rownum,]
  gene_names <- d1["gene_n"] %>%
    stringi::stri_split_regex(pattern = ",") %>%
    base::unlist()

  gfc_means <- gfc_dat[gfc_dat$Gene %in% gene_names,] %>%
    dplyr::select(-Gene) %>%
    base::colMeans()

  d1$conditions <- base::paste0(base::names(gfc_means), collapse = "#")
  d1$grp_means <- base::paste0(base::round(gfc_means,3) , collapse = ",")
  return(d1)

}


#' Function That Reads Gene Expression Matrices
#' 
#' Reads in the gene expression data from a matrix format. Files can be provided in most common file formats (.txt, .csv).
#'  Rows must correspond to genes, columns must correspond to samples.
#' @param file A string defining the file path.
#' @param rown A Boolean. Whether or not the file has rownames. Default is TRUE.
#' @param sep The separator of the file. Default is "\t" for tab separated files.
#' @param gene_symbol_col A String. Name of the column that contains the gene symbols.
#' @noRd

# function to read expression files:
read_expression_data <- function(file, rown = T, sep = "\t", gene_symbol_col){
  
  if(rown){
    
    expression_data <- utils::read.table(file = file, row.names = 1,
                                  stringsAsFactors = F, sep = sep, check.names = F, header = T, quote = "")
    
  }else{
    
    expression_data <- utils::read.table(file = file, 
                                  stringsAsFactors = F, sep = sep, check.names = F, header = T, quote = "")
  }
  expression_data <- make_rownames_unique(counts = expression_data, gene_symbol_col = gene_symbol_col)
  
  return(expression_data)
}


#' Remove Rows With Duplicate Rownames
#' 
#' The function removes all duplicate genes and only keeps the first occurance. It also drops all non-numeric columns.
#' @param counts The count matrix.
#' @param gene_symbol_col A String. Name of the column that contains the gene symbols.
#' @noRd

make_rownames_unique <- function(counts, gene_symbol_col){
  counts <- counts[!base::duplicated(counts[gene_symbol_col]),] %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., gene_symbol_col)
  
  # remove all non-numeric columns (description, gene id, etc.):
  for (x in base::colnames(counts)){
    
    if(!base::is.numeric(counts[[x]])){
      
      counts[[x]] <- NULL
      
    }
    
  }
  
  if(base::ncol(counts) == 0){
    
    print("All columns were deleted when removing non-numeric columns. Please check the data type of your expression values.")
    
  }
  
  return(counts)
}



#' Function That Reads The Annotation File
#' 
#' Reads in the annotation data from a matrix format. Files can be provided in most common file formats (.txt, .csv).
#'  Rows must correspond to samples, columns must correspond to meta information categories.
#'  Transforms all columns to factors.
#' @param file A string defining the file path.
#' @param rown A Boolean. Whether or not the file has rownames. Default is TRUE.
#' @param sep The separator of the file. Default is "\t" for tab separated files.
#' @param sample_col A String. Name of the column that contains the sample IDs.
#' @noRd

read_anno <- function(file, rown = T, sep = "\t", sample_col){
  
  if(rown){
    
    anno <- utils::read.table(file = file, row.names = 1,
                       stringsAsFactors = F, sep = sep, check.names = T, header = T) 
    
  }else{
    anno <- utils::read.table(file = file, 
                       stringsAsFactors = F, sep = sep, check.names = T, header = T)
    
  }
  base::rownames(anno) <- dplyr::pull(anno, sample_col) 
  
  anno[] <- base::lapply(anno, base::factor)
  
  return(anno)
}

#' Internal Implementation Of run_expression_analysis_1()
#' 
#' @param x An Integer. Gives the dataset that is currently processed.
#' @param padj A String. Defines the method to be used for p-value adjustment. Valid values are "none" (default) or values for "method" in stats::p.adjust.
#' @param export A Boolean. If TRUE, correlation values and p-values will be exported. 
#'  This can save time if you plan on re-running the analysis since computing pari-wise correlations is a bottleneck of the analysis. Default is FALSE.
#' @param import A list. Each slot represents a dataset and contains a vector of two file names, the first is the name of the correlation file of that dataset exported in previous runs.
#'  The second is the name of the p-value file exported in a previous run. For details, see the information file found in the repository. Default is NULL.
#' @param bayes Sánchez-Taltavull et al. (2016) suggest superiority of Bayesian correlation analysis to Pearson correlation in some cases. 
#'  Therefore, the Pearson correlation values can be weighted with Bayesian correlation values. To do so, set the “bayes”-parameter to TRUE. Default is FALSE, using only Pearson correlations.
#' @param alpha A numeric value from [0,1]. Allows to adjust the strength of the Bayes weighting: For alpha = 0 the Pearson correlation values remain unaltered, for alpha = 1 the Pearson correlation value and the Bayesian correlation value contribute equally to the final correlation.
#' @param prior An integer, either 2 or 3, using prior 2 or 3 for the Bayes weighting as described in "Bayesian correlation analysis for sequence count data" by Sánchez-Taltavull et al. (2016).
#' @param corr_method Method for the correlation calculation if bayes is FALSE. Either 'pearson', 'spearman' or 'rho'(single-cell).
#' @noRd


run_expression_analysis_1_body <- function(x, 
  bayes, 
  prior, 
  alpha, 
  padj, 
  export, 
  import,
  corr_method){

  message("Currently processed dataset: ", hcobject[["layers_names"]][x])
  output <- list()
  
  # retrieve gene count matrix for current layer:
  count_table <- hcobject[["data"]][[base::paste0("set", x, "_counts")]]

  # retrieve meta information for current layer:
  anno_table <- hcobject[["data"]][[base::paste0("set", x, "_anno")]]

  # filter for top most variant genes if required:

  if(!hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]] == "all"){
    message("...extracting ", hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]], " top variant genes...")
  }
  
    # sort the count data based on the genes' decreasing variance:
  ds = count_table[base::order(base::apply(count_table,1, stats::var), decreasing=T),]
  
  output[["ds"]] <- ds

    # filter:
  if(!hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]] == "all"){
    dd2 <- utils::head(ds, hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]])
  }else{
    dd2 <- ds
  }
  
  output[["topvar"]] <- dd2
  dd2 = t(dd2)
  
  
  # calculate pair-wise correlations:
  corr_calc_out <- pwcorr(dd2 = dd2, 
                layer_set = hcobject[["layer_settings"]][[base::paste0("set", x)]],
                bayes = bayes, 
                prior = prior, 
                alpha = alpha, 
                padj = padj, 
                export = export, 
                layer = x, 
                import = import,
                corr_method = corr_method)

  output[["corr_calc_out"]] <- corr_calc_out
  
  # calculate cut-off statistics:
  message("...calculating cutoff statistics...")
  cutoff_stats <- base::do.call("rbind", base::lapply(X = corr_calc_out[["range_cutoff"]],
                                         FUN = cutoff_prep,
                                         corrdf_r = corr_calc_out[["correlation_df_filt"]],
                                         print.all.plots = hcobject[["layer_settings"]][[base::paste0("set", x)]][["print_distribution_plots"]],
                                         x = x))
  
  
  output[["cutoff_stats"]] <- cutoff_stats
  
  # reshape cutoff stats:
  cutoff_calc_out <- reshape_cutoff_stats(cutoff_stats = cutoff_stats)
  
  output[["cutoff_calc_out"]] <- cutoff_calc_out
  knitr::kable(cutoff_calc_out[["cutoff_stats_concise"]], caption = "Correlation cut-off stats")  
  
  
  return(output)
  
}


#' Calculate p-value
#' 
#' Calculates the p-value ofd a given value based on a reference distribution.
#' @param x Value for which to calculate the p-value
#' @param mu The mean of the reference population
#' @param sigma The standard deviation of the reference population
#' @param n The size of the reference population
#' @noRd

calc_pval <- function(x, mu, sigma, n){
  z <- (x-mu)/(sigma/base::sqrt(n))
  p <- stats::pnorm(-base::abs(z))
  return(p)
}




#' Calculate Pair-Wise Correlations
#' 
#' The function calculates the pair-wise correlations for all genes in each dataset and performs multiple-testing correction. 
#'  Alternatively, previously calculated correlations and p-values can be imported. Pearson correlation coefficient may be weighted with Bayesian correlations.
#'  Correlations and p-values may be exported.
#' @param dd2 Transposed gene expression matrix.
#' @param padj A String. Defines the method to be used for p-value adjustment. Valid values are "none" (default) or values for "method" in stats::p.adjust.
#' @param export A Boolean. If TRUE, correlation values and p-values will be exported. 
#'  This can save time if you plan on re-running the analysis since computing pari-wise correlations is a bottleneck of the analysis. Default is FALSE.
#' @param import A list. Each slot represents a dataset and contains a vector of two file names, the first is the name of the correlation file of that dataset exported in previous runs.
#'  The second is the name of the p-value file exported in a previous run. For details, see the information file found in the repository. Default is NULL.
#' @param bayes Sánchez-Taltavull et al. (2016) suggest superiority of Bayesian correlation analysis to Pearson correlation in some cases. 
#'  Therefore, the Pearson correlation values can be weighted with Bayesian correlation values. To do so, set the “bayes”-parameter to TRUE. Default is FALSE, using only Pearson correlations.
#' @param alpha A numeric value from [0,1]. Allows to adjust the strength of the Bayes weighting: For alpha = 0 the Pearson correlation values remain unaltered, for alpha = 1 the Pearson correlation value and the Bayesian correlation value contribute equally to the final correlation.
#' @param prior An integer, either 2 or 3, using prior 2 or 3 for the Bayes weighting as described in "Bayesian correlation analysis for sequence count data" by Sánchez-Taltavull et al. (2016).
#' @param layer_set The layer specific settings for this layer.
#' @param layer An Integer indicating the currently processed dataset.
#' @noRd

pwcorr <- function(dd2, 
  layer_set, 
  bayes, 
  prior, 
  alpha, 
  padj, 
  export, 
  layer, 
  import,
  corr_method){ 
  
  
  message("...calculating pairwise correlations...")
  
  output <- list()
  
  # import of pre-calculated correlation values and their p-values:
  if(base::length(import) > 1){
    if(!base::is.na(import[layer])){
      # import matrix
      message("...importing correlation matrix from file...")
      correlation_matrix <- list()
      correlation_matrix[["r"]] <- readr::read_table2(import[[layer]][1]) %>% base::as.matrix()
      correlation_matrix[["P"]] <- readr::read_table2(import[[layer]][2]) %>% base::as.matrix()
      base::rownames(correlation_matrix[["r"]]) <- base::colnames(correlation_matrix[["r"]])
      base::rownames(correlation_matrix[["P"]]) <- base::colnames(correlation_matrix[["P"]])
    }else{
      if(corr_method == 'rho'){
        correlation_matrix <- list()
        correlation_matrix[["r"]] <- propr::perb(counts = base::as.matrix(dd2), select = colnames(dd2))@matrix
        pmat <- lapply(as.vector(correlation_matrix[["r"]]), 
                       calc_pval, 
                       mu = mean(correlation_matrix[["r"]]),
                       sigma = sd(correlation_matrix[["r"]]), 
                       n = ncol(correlation_matrix[["r"]])*nrow(correlation_matrix[["r"]])) %>% 
          unlist() %>% 
          matrix(., nrow = nrow(correlation_matrix[["r"]]), byrow = FALSE)
        colnames(pmat) <- colnames(correlation_matrix[["r"]])
        rownames(pmat) <- rownames(correlation_matrix[["r"]])
        correlation_matrix[["P"]] <- pmat
      }else{
        correlation_matrix <- Hmisc::rcorr(base::as.matrix(dd2), type=corr_method)
      }
    }
  }else if(base::length(import) == 1 & !base::is.null(import)){
    # import matrix
    message("...importing correlation matrix from file...")
    correlation_matrix <- list()
    correlation_matrix[["r"]] <- readr::read_table2(import[[layer]][1]) %>% base::as.matrix()
    correlation_matrix[["P"]] <- readr::read_table2(import[[layer]][2]) %>% base::as.matrix()
    base::rownames(correlation_matrix[["r"]]) <- base::colnames(correlation_matrix[["r"]])
    base::rownames(correlation_matrix[["P"]]) <- base::colnames(correlation_matrix[["P"]])
  }else{
    if(corr_method == 'rho'){
      correlation_matrix <- list()
      correlation_matrix[["r"]] <- propr::perb(counts = base::as.matrix(dd2), select = colnames(dd2))@matrix
      pmat <- lapply(as.vector(correlation_matrix[["r"]]), 
                     calc_pval, 
                     mu = mean(correlation_matrix[["r"]]),
                     sigma = sd(correlation_matrix[["r"]]), 
                     n = ncol(correlation_matrix[["r"]])*nrow(correlation_matrix[["r"]])) %>% 
        unlist() %>% 
        matrix(., nrow = nrow(correlation_matrix[["r"]]), byrow = FALSE)
      colnames(pmat) <- colnames(correlation_matrix[["r"]])
      rownames(pmat) <- rownames(correlation_matrix[["r"]])
      correlation_matrix[["P"]] <- pmat
    }else{
      correlation_matrix <- Hmisc::rcorr(base::as.matrix(dd2), type=corr_method)
    }
  }
  
  # export correlations and p-values for future re-runs to avoid the bottleneck:
  if(export){
    utils::write.table(x = correlation_matrix[["r"]], 
                file = base::paste0(hcobject[["working_director"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/correlation_matrix_", hcobject[["layers_names"]][layer], "_correlations.txt"), 
                quote = F, row.names = F, col.names = T, dec = ".")
    utils::write.table(x = correlation_matrix[["P"]], 
                file = base::paste0(hcobject[["working_director"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/correlation_matrix_", hcobject[["layers_names"]][layer], "_pvalues.txt"), 
                quote = F, row.names = F, col.names = T, dec = ".")
  }

  # reshaping the correlation amtrix to a dataframe with 4 columns (gene1, gene2, correaltion value, p-value):
  ind <- base::which(base::upper.tri(correlation_matrix[["r"]], diag = F), arr.ind = TRUE)
  correlation_df <- base::cbind(ind, correlation_matrix[["r"]][ind]) %>% base::as.data.frame()
  correlation_df[,1] <- base::colnames(dd2)[correlation_df[,1]]
  correlation_df[,2] <- base::colnames(dd2)[correlation_df[,2]]
  base::colnames(correlation_df) <- c("V1", "V2", "rval")
  correlation_df[["pval"]] <- correlation_matrix[["P"]][base::upper.tri(correlation_matrix[["P"]], diag = F)]
  
  # Bayes weighting:
  if(bayes){
    correlation_df[["rval"]] <- bayes_weighting(dd2 = dd2, alpha = alpha, prior = prior, pearson_rval = correlation_df[["rval"]])
  }
    
  # multiple testing correction:
  if(!padj == "none"){
    message("...conducting multiple testing correction using method ", padj, "...")
    correlation_df[["pval"]] <- stats::p.adjust(correlation_df[["pval"]], method = padj)
  }
  
  
  #retain rows which have pval (adj) < 0.05, and correlations above 0
  correlation_df_filt <- correlation_df[correlation_df[["pval"]] < 0.05 & correlation_df[["rval"]] > 0,]
  
  
  #range of cutoff min to max (correlation)
  range_cutoff <- base::seq(from = layer_set[["min_corr"]] , to = base::max(correlation_df[["rval"]]) , 
                    length.out = layer_set[["range_cutoff_length"]])
  range_cutoff <- base::round(range_cutoff, 3)
  if(base::length(range_cutoff) > layer_set[["range_cutoff_length"]]) {
    range_cutoff <- range_cutoff[1:layer_set[["range_cutoff_length"]]]
  } else {
    range_cutoff <- range_cutoff
  }
  
  
  output[["corr_mat"]] <- correlation_matrix
  output[["correlation_df"]] <- correlation_df
  output[["correlation_df_filt"]] <- correlation_df_filt
  output[["range_cutoff"]] <- range_cutoff
  
  return(output)
}


#' Function That Weights Given Pearson Correlation Coefficients With Bayes Correlation Values
#' 
#' @param dd2 Transposed gene expression matrix.
#' @param alpha A numeric value from [0,1]. Allows to adjust the strength of the Bayes weighting: For alpha = 0 the Pearson correlation values remain unaltered, for alpha = 1 the Pearson correlation value and the Bayesian correlation value contribute equally to the final correlation.
#' @param prior An integer, either 2 or 3, using prior 2 or 3 for the Bayes weighting as described in "Bayesian correlation analysis for sequence count data" by Sánchez-Taltavull et al. (2016).
#' @param pearson_rval A vector of Pearson Correlation Coefficients that are to be weighted.
#' @noRd

bayes_weighting <- function(dd2, 
              alpha, 
              prior, 
              pearson_rval){
  # use prior 2 as described in "Bayesian correlation analysis for sequence count data" by Sánchez-Taltavull et al. (2016).
  if(prior == 2){
    message("---weighting pearson with bayes using prior: 2---")
    # bayes:
    corr_b <- Bayes_Corr_Prior2(X = base::t(dd2))
  # use prior 3 as described in "Bayesian correlation analysis for sequence count data" by Sánchez-Taltavull et al. (2016).
  }else{
    message("---weighting pearson with bayes using prior: 3---")
    # bayes:
    corr_b <- Bayes_Corr_Prior3(X = base::t(dd2))
  }

  # reshape matrix:
  ind <- base::which(base::upper.tri(corr_b, diag = F), arr.ind = TRUE)
  correlation_df2 <- base::cbind(ind, corr_b[ind]) %>% base::as.data.frame()
  correlation_df2[,1] <- base::colnames(dd2)[correlation_df2[,1]]
  correlation_df2[,2] <- base::colnames(dd2)[correlation_df2[,2]]
  base::colnames(correlation_df2) <- c("V1", "V2", "rval")
  weighted_pearson <- (pearson_rval + (alpha * correlation_df2[["rval"]]))/(1 + alpha)

  return(weighted_pearson)

}

#' Computing The Bayesian Correlations Assuming Second (Dirichlet-marginalized) Prior
#' 
#' This function has been taken from  "Bayesian correlation analysis for sequence count data" by Sánchez-Taltavull et al. (2016), https://doi.org/10.1371/journal.pone.0163595.s001
#' @noRd

Bayes_Corr_Prior2 <- function(X){
  d <- base::dim(X)
  alpha0 <- base::rep(1/d[1],d[2])
  beta0 <- base::rep(1-1/d[1],d[2])
  Bcorrvals <- Bayes_Corr(alpha0,beta0,X)
  return(Bcorrvals)
}


#' Computing The Bayesian Correlations Assuming Third (zero count-motivated) Prior
#' 
#' This function has been taken from  "Bayesian correlation analysis for sequence count data" by Sánchez-Taltavull et al. (2016), https://doi.org/10.1371/journal.pone.0163595.s001
#' @noRd

Bayes_Corr_Prior3 <- function(X){
  d <- base::dim(X)
  cs <- base::colSums(X)
  alpha0 <- (cs+1)/(base::max(cs)+1)
  beta0 <- base::rep(1,d[2])
  Bcorrvals <- Bayes_Corr(alpha0,beta0,X)
  return(Bcorrvals)
}

#' Function To Compute Bayesian Correlation Coefficients
#' 
#' This function has been taken from  "Bayesian correlation analysis for sequence count data" by Sánchez-Taltavull et al. (2016), https://doi.org/10.1371/journal.pone.0163595.s001
#' @noRd

Bayes_Corr <- function(alpha0, beta0, X){
  nrowsX <- base::nrow(X)
  k <- base::ncol(X)
  cs <- base::colSums(X)
  alphas <- base::matrix(base::rep(alpha0,nrowsX), nrow=nrowsX, byrow=TRUE) + X
  betas  <- base::matrix(base::rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + base::matrix(base::rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  
  # First BIG product term for covariance formula
  Psi <- alphas/alphasPLUSbetas - base::matrix(base::rep(base::rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE) 
  
  # Covariance matrix
  cov_mtrx <- Psi %*% base::t(Psi) / k
  
  # Variances (this is a column vector of length = nrowsX)
  var_vec <- base::as.matrix( ( base::rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + base::rowSums(Psi^2) )/k )
  
  Bcorrvals <- cov_mtrx / base::sqrt( var_vec %*% base::t(var_vec) )
  base::diag(Bcorrvals) <- 1
  return(Bcorrvals)
}



#' Calculate Cutoff Statistics
#' 
#' @param cutoff The minimum correlation fro which to filter the data.
#' @param corrdf_r Data frame holding the correlation coefficients.
#' @param print.all.plots Boolean. Whether or not to print the degree distribution plots for all cutoffs.
#' @param x An Integer. Gives the dataset that is currently processed.
#' @noRd

cutoff_prep <- function(cutoff, corrdf_r, print.all.plots, x, min_nodes = hcobject[["global_settings"]][["min_nodes_number_for_network"]]){
  
  
  ###filter correlations above the cutoff
  filteredmatrix = corrdf_r[corrdf_r[["rval"]] >= cutoff,]
  
  ##create expected df with initial values
  output <- base:: data.frame(R.squared = 0,
                              degree = 0,
                              Probs = 0,
                              cutoff = cutoff,
                              no_edges = 0,
                              no_nodes = 0,
                              no_of_networks = 0)
  

  rownums <- base::ifelse(base::nrow(filteredmatrix)==0, "zero_rows", "gtzero")
  
  base::switch(rownums, zero_rows = {output},
         gtzero={rsquaredfun(graph_df = filteredmatrix, cutoff = cutoff, print.all.plots = print.all.plots, min_nodes = min_nodes)})
  
}


#' Function For Calculating Network Statistics
#' 
#' @param graph_df Description of the network in matrix format with two columns giving edges between nodes and a third column giving edge weights.
#' @param cutoff The minimum correlation for which the edges are filtered.
#' @noRd


rsquaredfun <- function(graph_df, cutoff, print.all.plots, min_nodes = hcobject[["global_settings"]][["min_nodes_number_for_network"]]){

  #igraph object
  filt_igraph = igraph::graph_from_data_frame(graph_df, directed = F, vertices = NULL)

  #number of components
  graph_components <- igraph::components(filt_igraph)
  num_networks = graph_components[["csize"]][graph_components[["csize"]] >= min_nodes] %>% 
    base::length()

  # remove nodes from too small components:
  gene_to_comp <- base::data.frame(gene = base::names(graph_components[["membership"]]), component = graph_components[["membership"]])
  comps_to_keep <- base::which(graph_components[["csize"]] >= min_nodes)
  nodes_to_remove <- dplyr::filter(gene_to_comp, !component %in% comps_to_keep) %>% dplyr::pull(., "gene")
  filt_igraph <- igraph::delete.vertices(filt_igraph, nodes_to_remove)

  #number of nodes
  num_nodes <- igraph::vcount(filt_igraph)
  #number of edges
  num_edges <- igraph::gsize(filt_igraph)

  
  ##calculate stats
  if(igraph::vcount(filt_igraph) == 0){
    R.square <- NA
    degree <- NA
    probability <- NA
  }else{
    d <- igraph::degree(filt_igraph, mode = "all")
    dd <- igraph::degree.distribution(filt_igraph, mode = "all", cumulative = FALSE)
    degree <- 1:base::max(d)
    probability <- dd[-1]
    nonzero.position <- base::which(probability != 0)
    probability <- probability[nonzero.position]
    degree <- degree[nonzero.position]
    
    if(base::length(probability) == 0){
      R.square <- 0
    }else{
      forplot <- base::data.frame(probability = probability, degree=degree)
      reg <- stats::lm(base::log(probability) ~ base::log(degree))
      cozf <- stats::coef(reg)
      power.law.fit <- function(x) base::exp(cozf[[1]] + cozf[[2]] * base::log(x))
      alpha <- -cozf[[2]]
      R.square <- base::summary(reg)$r.squared
      
      
      if(print.all.plots) {
        
        degree_distribution_wd <- base::paste0("dir_DegreeDistribution_", hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]])
        
        if(!degree_distribution_wd %in% base::list.dirs(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]]))) {
          base::dir.create(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], degree_distribution_wd))}
        
        dd_plot = ggplot2::ggplot(forplot, ggplot2::aes(x = base::log(degree), y = base::log(probability))) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method="lm") +
          ggplot2::theme_bw()
        
        ggplot2::ggsave(path = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], degree_distribution_wd),
               filename = base::paste0("Degree_distribution_plot_", cutoff,"_set_", x, ".pdf"),
               dd_plot, device = cairo_pdf)
        
        
      }
      
      
    }
  }
  
  output <- base::data.frame(R.squared= R.square,
                             degree = degree,
                             Probs = probability,
                             cutoff = cutoff,
                             no_edges = num_edges,
                             no_nodes = num_nodes,
                             no_of_networks = num_networks)
  
  return(output)
  
}

####### ADD FUNCTIONS FROM optimal_cutoff_MultiOmics.R #######
#' Reshape Cutoff Statistics
#' 
#' The function reshapes the cutoff statistics for downstream usage.
#' @param cutoff_stats The native representation of the cutoff statistics.
#' @noRd

reshape_cutoff_stats <- function(cutoff_stats){
  
  
  output <- list()
  
  cutoff_stats_concise <-  cutoff_stats %>% 
    dplyr::select(R.squared, cutoff, no_edges, no_nodes, no_of_networks) %>% 
    dplyr::distinct()
  
  cutoff_stats_concise <- cutoff_stats_concise %>% dplyr::filter(no_of_networks != 0)
  base::rownames(cutoff_stats_concise) <- cutoff_stats_concise[["cutoff"]]
  cutoff_stats_concise <- cutoff_stats_concise[,-2]
  
  
  
  crit_minmax = c("max","max", "max", "min" )
  base::names(crit_minmax) = base::colnames(cutoff_stats_concise)
  
  
  
  normalizationTypes <- base::rep("percentageOfMax", base::ncol(cutoff_stats_concise))
  base::names(normalizationTypes) = base::colnames(cutoff_stats_concise)
  if(base::nrow(cutoff_stats_concise) == 0){
    print("cutoff_stats_concise is empty")
    return(base::data.frame())
  }
  nPT <- MCDA::normalizePerformanceTable(cutoff_stats_concise[, c("R.squared", "no_edges", "no_nodes", "no_of_networks")], normalizationTypes)
  w <- base::c(0.5, 0.1, 0.5, -1)
  base::names(w) <- base::colnames(nPT)
  ws <- MCDA::weightedSum(nPT,w)
  ranked_ws <- base::rank(-ws) %>% base::sort()
  
  
  calculated_optimal_cutoff <- base::as.numeric(base::names(ranked_ws[1]))
  
  
  
  stats_calculated_optimal_cutoff <- cutoff_stats[cutoff_stats[["cutoff"]] == calculated_optimal_cutoff, c("degree", "Probs")]
  
  dd_plot_calculated_optimal <- ggplot2::ggplot(stats_calculated_optimal_cutoff, ggplot2::aes(x = base::log(degree), y = base::log(Probs))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method="lm") +
    ggplot2::theme_bw() + 
    ggplot2::ggtitle(base::paste0("Calculated optimal correlation cut-off [", calculated_optimal_cutoff, "]"))
  
 
  output[["cutoff_stats_concise"]] <- cutoff_stats_concise
  output[["dd_plot_calculated_optimal"]] <- dd_plot_calculated_optimal
  output[["optimal_cutoff"]] <- calculated_optimal_cutoff
  return(output)

}

#' Internal Realization Of run_expression_analysis_2()
#' 
#' Iterates over all datasets to perform the actions described in run_expression_analysis_2().
#' @param x An integer giving the number of the dataset currently processed.
#' @param grouping_v A string giving a column name present in all annotation files, if this variable shall be used for grouping the samles isntead of the variable of interest. Default is NULL.
#' @param plot_HM A Boolean. Whether or not to plot the heatmap (for networks with many genes this may be very demanding for your computer if you are running the analysis locally). Default is TRUE.
#' @param method The method used for clustering the heatmap in the pheatmap function. Default is "complete".
#' @param additional_anno A list, with one slot per data set. A slot contains a vector of column names from that data set’s annotation file that you wish to annotate with. 
#'  If for some of the data sets you don’t wish any further annotation, you can set the corresponding list slot to NULL. Default is NULL.
#' @noRd

run_expression_analysis_2_body <- function(x, grouping_v, plot_HM, method, additional_anno){
  
  message("...Currently processed dataset: ", hcobject[["layers_names"]][x], '...')

  hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part2"]][["heatmap_out"]] <<- heatmap_network_genes(x = x, 
                                                              plot_HM = plot_HM,
                                                              method = method, 
                                                              additional_anno = additional_anno)

  
  hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part2"]][["GFC_all_genes"]] <<- GFC_calculation(info_dataset = hcobject[["data"]][[base::paste0("set", x, "_anno")]],
                                                                           grouping_v = grouping_v,
                                                                           x = x)
  
}



#' Plot Heatmap Of Network Genes
#' 
#' Plots a heatmap of the genes present in the network, thus those genes left after filtering for the minimum correlation and removing too smal graph components.
#' @param x An integer giving the number of the dataset currently processed.
#' @param plot_HM A Boolean. Whether or not to plot the heatmap. 
#' @param method The method used for clustering in the pheatmap function.
#' @param additional_anno A vector of strings giving other columns from the annotation to be annotated in the heatmap in addition to the variable of interest.
#' @noRd

heatmap_network_genes <- function(x, plot_HM, method, additional_anno){


  # extract annotation data for current data layer:
  info_dataset <- hcobject[["data"]][[base::paste0("set",x,"_anno")]]

  message("...creating heatmap of network genes...")

  # collect function output:
  output <- list()

  # conditions to annotate in the heatmap:
  all_conditions <- base::unique(base::c(additional_anno, hcobject[["global_settings"]][["voi"]]))
  all_conditions <- all_conditions[!base::is.null(all_conditions)]
  
  # extract data exceeding set cutoff:
  filt_cutoff_data <- hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["corr_calc_out"]][["correlation_df_filt"]] %>% 
              dplyr::filter(., rval >= hcobject[["cutoff_vec"]][x])

  # build network based on filtered data:
  filt_cutoff_graph <- igraph::graph_from_data_frame(filt_cutoff_data, directed = FALSE)

  # remove too small components:
  graph_components <- igraph::components(filt_cutoff_graph)

  # remove nodes from too small components from the network:
  gene_to_comp <- base::data.frame(gene = base::names(graph_components[["membership"]]), component = graph_components[["membership"]])
  comps_to_keep <- base::which(graph_components[["csize"]] >= hcobject[["global_settings"]][["min_nodes_number_for_network"]])
  nodes_to_remove <- dplyr::filter(gene_to_comp, !component %in% comps_to_keep) %>% dplyr::pull(., "gene")
  filt_cutoff_graph <- igraph::delete.vertices(filt_cutoff_graph, nodes_to_remove)
  
  # remove edges from removed nodes from the cutoff data:
  filt_cutoff_data <- dplyr::filter(filt_cutoff_data, V1 %in% igraph::V(filt_cutoff_graph)$name & V2 %in% igraph::V(filt_cutoff_graph)$name)
  
  filt_cutoff_counts <- hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["ds"]][base::row.names(hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["ds"]]) %in% base::names(igraph::V(filt_cutoff_graph)),]
  corresp_info = info_dataset[base::rownames(base::t(hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["topvar"]])) %in% base::rownames(info_dataset),]
  
  output[["filt_cutoff_graph"]] <- filt_cutoff_graph
  output[["filt_cutoff_data"]] <- filt_cutoff_data
  
  
  message("...After using the optimal cutoff of ", hcobject[["cutoff_vec"]][x], " the number of edges = ", 
              base::nrow(filt_cutoff_data), " and the number of nodes = ", base::nrow(filt_cutoff_counts), '...')
 
  col_list <- list()
  for(i in all_conditions){
    tmp_col <- ggsci::pal_nejm(alpha = 1)(8)
    # if there are more groups than colours, expand palette:
    if(base::length(base::unique(hcobject[["data"]][[paste0("set", x, "_anno")]][,i])) > base::length(tmp_col)){
      tmp_col <- grDevices::colorRampPalette(tmp_col)(base::length(base::unique(hcobject[["data"]][[base::paste0("set", x, "_anno")]][,i])))
    }else{
      tmp_col <- tmp_col[1:base::length(base::unique(hcobject[["data"]][[base::paste0("set", x, "_anno")]][,i]))]
    }
    base::names(tmp_col) <- base::unique(hcobject[["data"]][[base::paste0("set", x, "_anno")]][,i])
    col_list[[i]] <- tmp_col
  }

  # filter annotation for pheatmap:
  anno_df <- dplyr::select(hcobject[["data"]][[base::paste0("set", x, "_anno")]], tidyselect::all_of(all_conditions))
  
  
  if(plot_HM){
    heatmap_filtered_counts <- pheatmap::pheatmap(mat = filt_cutoff_counts ,
                                                  color = base::rev(RColorBrewer::brewer.pal(11, "RdBu")),
                                                  scale = "row",
                                                  cluster_rows = T,
                                                  cluster_cols = T,
                                                  annotation_colors = col_list,
                                                  annotation_col = anno_df,
                                                  fontsize = 8,
                                                  show_rownames = F, 
                                                  show_colnames = T, 
                                                  annotation_names_col = T, 
                                                  clustering_distance_cols = "euclidean", 
                                                  clustering_method = method)
  }else{
    
    grDevices::pdf(file = NULL)
    heatmap_filtered_counts <- pheatmap::pheatmap(mat = filt_cutoff_counts ,
                                                  color = base::rev(RColorBrewer::brewer.pal(11, "RdBu")),
                                                  scale = "row",
                                                  cluster_rows = T,
                                                  cluster_cols = T,
                                                  annotation_colors = col_list,
                                                  annotation_col = anno_df,
                                                  fontsize = 8,
                                                  show_rownames = F, 
                                                  show_colnames = T, 
                                                  annotation_names_col = T, 
                                                  clustering_distance_cols = "euclidean", 
                                                  clustering_method = method)
    grDevices::dev.off()
  }
  
  output[["heatmap"]] <- heatmap_filtered_counts
  ggplot2::ggsave(filename = base::paste0("Heatmap_topvar_genes", x,".pdf"), plot = heatmap_filtered_counts, device = cairo_pdf,
         path = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]]), width = 7, height = 10, units = "in", limitsize = F)

  return(output)
}




#' Truncated Fold Changes
#'
#' Function calculates fold changes between given columns and a reference and truncates them if the FCs exceed the set maximum GFC range.
#' @param grp A vector of column names for which the fold changes shall be calculated.
#' @param group_means Name of the column containing the reference level from which the FC shall be calculated.
#' @param trans_norm The dataframe containing the data.
#' @noRd

gfc_calc <- function(grp, trans_norm, group_means){

  df1 <- trans_norm[,grp]
  df2 <- gtools::foldchange(df1, group_means) %>% base::ifelse(. >  hcobject[["global_settings"]][["range_GFC"]], hcobject[["global_settings"]][["range_GFC"]], . ) %>% 
    base::ifelse(. < (-hcobject[["global_settings"]][["range_GFC"]]), -hcobject[["global_settings"]][["range_GFC"]],.) %>%
    base::as.data.frame()
  base::colnames(df2) <-  base::paste0("", grp)
  return(df2)
}


#' Calculate Group Fold Changes
#' 
#' Function that computes fold changes for the define sample groups either with reference to a control group or to the mean of all groups.
#' @param info_dataset The annotation dataframe for the dataset that contains the group labels.
#' @param grouping_v A string giving a column name present in all annotation files, if this variable shall be used for grouping the samles isntead of the variable of interest.
#' @param x An integer giving the number of the curently processed dataset.
#' @noRd


GFC_calculation <- function(info_dataset, grouping_v, x) {
  
  message("...calculate Group-Fold-Changes...")
  
  if(!base::is.null(grouping_v)){
    info_dataset[["grpvar"]] <- info_dataset[,base::c(grouping_v)]
    print(base::paste0("User-defined variable ", grouping_v, " will be used for grouping the data."))
  }
  else if(base::intersect(hcobject[["global_settings"]][["voi"]], base::colnames(info_dataset)) %>% base::length() > 0){
    
    message("...Variable: '", hcobject[["global_settings"]][["voi"]],
                "' will be used as grouping variables...")
    
    info_dataset[["grpvar"]] <- purrr::pmap(info_dataset[base::intersect(hcobject[["global_settings"]][["voi"]], base::colnames(info_dataset))],
                                                        paste, sep = "-") %>% base::unlist()
  } else {
    
    message("...The first column in the metadata will be used as the grouping variable",
                "since the voi_id is not present in the metadata...")

    info_dataset[["grpvar"]] = info_dataset[,1]
    
  }
  
  
  if(hcobject[["global_settings"]][["control"]] == "none"){

      message("...GFC calculation with foldchange from mean...")

      # GFC calculation with foldchange from mean
      norm_data_anno <- base::merge(base::t(hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["topvar"]]), base::subset(info_dataset, select = base::c("grpvar")), by = "row.names", all.x = T)

      norm_data_anno <- norm_data_anno[,-1]
      
      norm_data_anno <- norm_data_anno[ , base::c(base::ncol(norm_data_anno) , 1:(base::ncol(norm_data_anno)-1))]
      
      trans_norm <- stats::setNames(base::data.frame(base::t(norm_data_anno[ , -1])) , norm_data_anno[,1])
      
      if(hcobject[["global_settings"]][["data_in_log"]] == TRUE){
        trans_norm <- antilog(trans_norm , 2)
      }
      
     
      
      trans_norm <- base::t(base::apply(trans_norm , 1 , function(x) base::tapply(x , base::colnames(trans_norm) , base::mean)))
      trans_norm <- base::cbind(trans_norm , base::rowMeans(trans_norm))
      
      base::colnames(trans_norm)[base::ncol(trans_norm)] <- "group_mean"
      grplist <- base::colnames(trans_norm)[-(base::ncol(trans_norm))]
      
      
      GFC_all_genes <- base::do.call("cbind", base::lapply(grplist, gfc_calc, trans_norm = trans_norm, group_means = trans_norm[,"group_mean"]))
      GFC_all_genes <- base::round(GFC_all_genes,3)
      GFC_all_genes$Gene <- base::rownames(GFC_all_genes)
      
      
      return(GFC_all_genes)

  }else{
    message("...GFC calculation with foldchange from control...")
    # GFC calculation with foldchange from control
    norm_data_anno = base::merge(base::t(hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["topvar"]]), base::subset(info_dataset, select = base::c("grpvar")), by = "row.names", all.x = T)
    
    
    #contains a column called Row.names (1st col)
    norm_data_anno = norm_data_anno[,-1]

    norm_data_anno <- norm_data_anno[ , base::c(base::ncol(norm_data_anno) , 1:(base::ncol(norm_data_anno)-1))]
    
    
    trans_norm <- stats::setNames(base::data.frame(base::t(norm_data_anno[ , -1])) , norm_data_anno[,1])
    
    if(hcobject[["global_settings"]][["data_in_log"]] == TRUE){
      trans_norm[trans_norm == 0] <- 1
      trans_norm <- antilog(trans_norm , 2)
    }
    
    trans_norm <- base::t(base::apply(trans_norm , 1 , function(x) base::tapply(x , base::colnames(trans_norm) , base::mean)))
    
    
    
    
    trans_norm_no_ctrl <- trans_norm[, base::grepl(hcobject[["global_settings"]][["control"]], base::colnames(trans_norm), ignore.case = T) == F]
    if(base::is.vector(trans_norm_no_ctrl)){
      
      trans_norm_no_ctrl <- base::data.frame(trans_norm_no_ctrl = trans_norm_no_ctrl)
      
      base::colnames(trans_norm_no_ctrl) <- base::colnames(trans_norm)[base::grepl(hcobject[["global_settings"]][["control"]], 
                                                                 base::colnames(trans_norm), ignore.case = T) == F]
      base::rownames(trans_norm_no_ctrl) <- base::rownames(trans_norm)
    }

    trans_norm <- base::cbind(trans_norm_no_ctrl, trans_norm[, base::grepl(hcobject[["global_settings"]][["control"]], 
                                                               base::colnames(trans_norm), ignore.case = T) == T])
    
    base::rownames(trans_norm) <- base::rownames(trans_norm_no_ctrl)
    base::colnames(trans_norm)[base::ncol(trans_norm)] <- "group_mean"
    grplist <- base::colnames(trans_norm)[-(base::ncol(trans_norm))]
    
    GFC_all_genes <- base::do.call("cbind", base::lapply(grplist, gfc_calc, trans_norm = trans_norm, group_means = trans_norm[,"group_mean"]))
    GFC_all_genes <- base::round(GFC_all_genes, 3)
    base::rownames(GFC_all_genes) <- base::rownames(trans_norm)
    GFC_all_genes$Gene = base::rownames(GFC_all_genes)
    tmp_col_names <- base::colnames(GFC_all_genes)[base::grepl(hcobject[["global_settings"]][["control"]], 
                                                   base::colnames(GFC_all_genes), ignore.case = T)== F]
    GFC_all_genes <- GFC_all_genes[, base::colnames(GFC_all_genes) %in% tmp_col_names]
    
    return(GFC_all_genes)
  }
  
  
}



#' Get Union Of Graphs
#' 
#' The function builds a multigraph from the union of all layer-specific networks. 
#'  Thus, also network parts that are unique to some layers will be present in the resulting integrated network, creating an integrated network that provides a wholistic, cross-dataset view on gene co-expression.
#' @noRd

get_union <- function(){
  
  combined_edgelist <- NULL
  
  for (x in 1:base::length(hcobject[["layer_specific_outputs"]])){
    combined_edgelist <- base::rbind(combined_edgelist, hcobject[["layer_specific_outputs"]][[x]][["part2"]][["heatmap_out"]][["filt_cutoff_data"]])
  }
  combined_edgelist$weight <- combined_edgelist$rval
  combined_edgelist$rval <- NULL
  combined_edgelist$pval <- NULL
  

  hcobject[["integrated_output"]][["combined_edgelist"]] <<- combined_edgelist
}


#' Get Intersection With Reference Graph
#' 
#' It generates a multigraph of the reference network. 
#'  The vertices of the resulting network will be identical to the reference network, but the edges connecting them will be greatly impacted by the other datasets.
#' @param with Either an integer giving the number of the dataset to be used as reference (e.g., 1) or the name given to the layer.
#' @noRd

get_intersection <- function(with){

  # change to numeric representation of the dataset in case it was given as the name of a dataset:
  if(with %in% hcobject[["layers_names"]]){
    with <- match(with, hcobject[["layers_names"]])
  }

  # change to numeric representation of the dataset in case it was given as a string:
  if(startsWith(with, "set")){
    with <- base::strsplit(base::as.character(with), split = "set")[[1]][2] %>% base::as.numeric()
  }
  
  # get edgelist of the reference network:
  combined_edgelist <- hcobject[["layer_specific_outputs"]][[with]][["part2"]][["heatmap_out"]][["filt_cutoff_data"]]

  combined_edgelist$merged <- base::paste0(combined_edgelist$V1 %>% base::as.character(), 
                                     combined_edgelist$V2 %>% base::as.character())
  combined_edgelist$revmerged <- base::paste0(combined_edgelist$V2 %>% base::as.character(), 
                                     combined_edgelist$V1 %>% base::as.character())
  # iterate over datasets:
  for (x in 1:base::length(hcobject[["layer_specific_outputs"]])){

    if(!x == with){
      # get edgelist of current dataset:
      tmp <- hcobject[["layer_specific_outputs"]][[x]][["part2"]][["heatmap_out"]][["filt_cutoff_data"]]

      tmp$merged <- base::paste0(tmp$V1 %>% base::as.character(), 
                           tmp$V2 %>% base::as.character())

      # check which edges overlap with reference network:
      tmp <- dplyr::filter(tmp, merged %in% combined_edgelist$merged | merged %in% combined_edgelist$revmerged)
      combined_edgelist <- base::rbind(combined_edgelist, tmp)
    }
  }
  combined_edgelist$weight <- combined_edgelist$rval
  combined_edgelist$rval <- NULL
  combined_edgelist$pval <- NULL
  combined_edgelist$merged <- NULL
  combined_edgelist$revmerged <- NULL
  
  hcobject[["integrated_output"]][["combined_edgelist"]] <<- combined_edgelist
}

#' Merge GFCs From Different Datasets
#' 
#' Combines the GFC values per gene from the different datasets and substitutes missing values.
#' @param GFC_when_missing The GFC value to enter when a gene has not been measured in a dataset, but in others. Default is the lower bound of the set GFC range.
#' @noRd

merge_GFCs <- function(GFC_when_missing = -hcobject[["global_settings"]][["range_GFC"]]){


  col_names_new_GFC <- NULL
  for (x in 1:base::length(hcobject[["layer_specific_outputs"]])){
    new_col <- base::colnames(hcobject[["layer_specific_outputs"]][[x]][["part2"]][["GFC_all_genes"]])[-base::ncol(hcobject[["layer_specific_outputs"]][[x]][["part2"]][["GFC_all_genes"]])]

    col_names_new_GFC <- base::c(col_names_new_GFC, new_col)
  }
  new_GFC <- NULL
  for(y in igraph::get.vertex.attribute(hcobject[["integrated_output"]][["merged_net"]])$name){
    
    line <- NULL
    
    
    for(z in 1:base::length(hcobject[["layer_specific_outputs"]])){
      if(y %in% hcobject[["layer_specific_outputs"]][[z]][["part2"]][["GFC_all_genes"]][["Gene"]]){

        GFC_tmp <- hcobject[["layer_specific_outputs"]][[z]][["part2"]][["GFC_all_genes"]]
        
        line <- base::c(line, GFC_tmp[GFC_tmp$Gene == y, base::colnames(GFC_tmp)[1:(base::ncol(GFC_tmp)-1)]]) %>%
          base::unlist(.)
        
      }else{
        line <- base::c(line, base::rep(GFC_when_missing, (base::length(base::colnames(hcobject[["layer_specific_outputs"]][[z]][["part2"]][["GFC_all_genes"]]))-1)))
      }
    }
    
    line <- base::as.data.frame(line) %>% base::t()
    new_GFC <- base::rbind(new_GFC, line)
  }
  new_GFC <- base::as.data.frame(new_GFC)
  
  base::colnames(new_GFC) <- base::c(col_names_new_GFC)
  new_GFC$Gene <- igraph::get.vertex.attribute(hcobject[["integrated_output"]][["merged_net"]])$name
  rownames(new_GFC) <- new_GFC$Gene
  return(new_GFC)
}


#' Fixes Missing Slashes
#' 
#' Adds slash to the end of provided directory if it is missing and gives and error if the provided directory does not exist.
#' @param directory A string describing a directory path.
#' @noRd

fix_dir <- function(directory){
  if(directory == FALSE){
    return(directory)
  }
  # check existance of path:
  if(!base::dir.exists(directory)){
    stop( "The directory '", directory, "' does not exist.")
  }
  # add slash if missing at the end:
  split_dir <- base::strsplit(directory, split = "")[[1]]
  if(split_dir[base::length(split_dir)] == "/"){
    return(directory)
  }else{
    return(base::paste0(directory, "/"))
  }
}


#' Non-Interactive Plot of Cutoff Statistics
#' 
#' Plots the R-squared vlaue, number of edges, number of genes and number of networks for different cut-offs as a ggplot.
#' @param cutoff_stats A dataframe of cutoff statistics generated in previous steps.
#' @param hline A list with four slots ("R.squared", "no_edges", "no_nodes", "no_networks") each of which can be set either to NULL (default) or a number to introduce a horizontal line for orientation at that value in the respective plot.
#' @param x An integer giving the number of the currently processed layer.
#' @noRd

plot_cutoffs_internal_static <- function(cutoff_stats,
                                  hline = list("R.squared" = NULL, "no_edges" = NULL, "no_nodes" = NULL, "no_networks" = NULL),
                                  x){
  cutoff_stats[["corr"]] <- base::rownames(cutoff_stats) %>% base::as.numeric()

  # plot R-squared values:
  p1 <- ggplot2::ggplot(cutoff_stats, ggplot2::aes(x = corr))+
    ggplot2::geom_vline(xintercept = base::c(base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.01)), color = "darkgrey")+
    ggplot2::geom_hline(yintercept = hline[[1]], color = "darkgrey", size = 1.7)+
    ggplot2::geom_line(ggplot2::aes(y = R.squared), color = "#374E55FF")+
    ggplot2::geom_point(ggplot2::aes(y = R.squared), color = "#374E55FF", size = 2)+
    ggplot2::theme_light()+
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(breaks = base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.005))+
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::ggtitle(hcobject[["layers_names"]][x])

  # plot number of edges:
  p2 <- ggplot2::ggplot(cutoff_stats, ggplot2::aes(x = corr))+
    ggplot2::geom_vline(xintercept = base::c(base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.01)), color = "darkgrey")+
    ggplot2::geom_hline(yintercept = hline[[2]], color = "darkgrey", size = 1.7)+
    ggplot2::geom_line(ggplot2::aes(y = no_edges), color = "#DF8F44FF")+
    ggplot2::geom_point(ggplot2::aes(y = no_edges), color = "#DF8F44FF", size = 2)+
    ggplot2::theme_light()+
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank())+
    ggplot2::scale_x_continuous(breaks = base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.005))+
    ggplot2::scale_y_continuous(labels = scales::comma)

  # plot number of nodes:
  p3 <- ggplot2::ggplot(cutoff_stats, ggplot2::aes(x = corr))+
    ggplot2::geom_vline(xintercept = base::c(base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.01)), color = "darkgrey")+
    ggplot2::geom_hline(yintercept = hline[[3]], color = "darkgrey", size = 1.7)+
    ggplot2::geom_line(ggplot2::aes(y = no_nodes), color = "#00A1D5FF")+
    ggplot2::geom_point(ggplot2::aes(y = no_nodes), color = "#00A1D5FF", size = 2)+
    ggplot2::theme_light()+
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank())+
    ggplot2::scale_x_continuous(breaks = base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.005))+
    ggplot2::scale_y_continuous(labels = scales::comma)

  # plot number of networks:
  p4 <- ggplot2::ggplot(cutoff_stats, ggplot2::aes(x = corr))+
    ggplot2::geom_vline(xintercept = base::c(base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.01)), color = "darkgrey")+
    ggplot2::geom_hline(yintercept = hline[[4]], color = "darkgrey", size = 1.7)+
    ggplot2::geom_point(ggplot2::aes(y = no_of_networks), color = "#B24745FF", size = 2)+
    ggplot2::scale_x_continuous(breaks = base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.005))+
    ggplot2::theme_light()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90))


  p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1, align = "v")
  return(p)
}

#' Interactive Plot of Cutoff Statistics
#' 
#' Plots the R-squared vlaue, number of edges, number of genes and number of networks for different cut-offs as an interactive widget using plotly.
#' @param cutoff_stats A dataframe of cutoff statistics generated in previous steps.
#' @param x An integer giving the number of the currently processed layer.
#' @noRd

plot_cutoffs_internal_interactive <- function(cutoff_stats, 
                            x){

  cutoff_stats[["corr"]] <- base::rownames(cutoff_stats) %>% base::as.numeric()
  
  # plot R-squared value:
  p1 <- plotly::plot_ly(cutoff_stats, x = ~corr, y = ~R.squared, type = 'scatter', 
                mode = 'lines+markers', name = "R²", line = list(color = "#374E55FF"), marker = list(color = "#374E55FF")) 
  # plot number of edges:
  p2 <- plotly::plot_ly(cutoff_stats, x = ~corr, y = ~no_edges, type = 'scatter', 
                mode = 'lines+markers', name = "no. edges", line = list(color = "#DF8F44FF"), marker = list(color = "#DF8F44FF"))
  # plot number of nodes:
  p3 <- plotly::plot_ly(cutoff_stats, x = ~corr, y = ~no_nodes, type = 'scatter', 
                mode = 'lines+markers', name = "no. nodes", line = list(color = "#00A1D5FF"), marker = list(color = "#00A1D5FF"))
  # plot number of networks:
  p4 <- plotly::plot_ly(cutoff_stats, x = ~corr, y = ~no_of_networks, type = 'scatter', 
                mode = "markers", name = "no. networks", marker = list(color = "#B24745FF"))
  # combine plots:
  p <- plotly::subplot(p1, p2, p3, p4, nrows = 4, shareX = T)
  
  # adjust layout:
  steps <- list()
  for(i in 1:base::length(cutoff_stats[["corr"]])){
    
    step <- list(args = list("marker.color",list(base::rep("#374E55FF", base::length(cutoff_stats[["corr"]])),
                                                 base::rep("#DF8F44FF", base::length(cutoff_stats[["corr"]])),
                                                 base::rep("#00A1D5FF", base::length(cutoff_stats[["corr"]])),
                                                 base::rep("#B24745FF", base::length(cutoff_stats[["corr"]])))), 
                 label = base::paste0(base::as.character(cutoff_stats[["corr"]][i]), ", R²: ", 
                                base::round(cutoff_stats[["R.squared"]][i], 3),
                                "; no. edges: ", cutoff_stats[["no_edges"]][i], "; no. nodes: ", 
                                cutoff_stats[["no_nodes"]][i], "; no. networks: ", cutoff_stats[["no_of_networks"]][i]), 
                 method = "restyle"
                 
    )
    
    for(j in 1:4){
      step[["args"]][[2]][[j]][i] <- "red"
    }
    
    steps[[i]] <- step
  }
  
  
  p <- p %>% plotly::layout(hovermode = "x unified") %>%
    plotly::layout(title = base::paste0("Cut-off selection guide: ", hcobject[["layers_names"]][x]),
           sliders = list(
             list(pad = list(t=60),
                  active = 2, 
                  currentvalue = list(prefix = "Cut-off: ", font = list(color = "black", size = 14)), 
                  steps = steps,
                  font = list(color = "white", size = 0))),
           shapes = list(line))
  return(p)
}

#' Helper Function To Plot TF Enrichment
#' @noRd

plot_TF <- function(hubs_df){
  #get column order from module heatmap:
  co <- ComplexHeatmap::column_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]]) %>% base::unlist() %>% base::as.numeric()
  co <- base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]][,-base::ncol(hcobject[["integrated_output"]][["GFC_all_layers"]])])[co]
  
  final_colnames <- NULL
  hub_exp <- base::apply(hubs_df, 2, function(x){
    genes <- x[!x == " "]
    if(base::length(genes) == 0){
      return(NULL)
    }
    
    out <- base::lapply(1:base::length(hcobject[["layers"]]), function(l){
      exp <- hcobject[["data"]][[base::paste0("set", l, "_counts")]][genes,]
      mean_exp <- base::lapply(base::unique(hcobject[["data"]][[base::paste0("set", l, "_anno")]][[hcobject[["global_settings"]][["voi"]]]]), function(c){
        samples <- base::rownames(hcobject[["data"]][[base::paste0("set", l, "_anno")]][hcobject[["data"]][[base::paste0("set", l, "_anno")]][[hcobject[["global_settings"]][["voi"]]]]==c,] )
        tmp <- dplyr::select(exp, samples) %>% base::rowMeans()
      }) %>% rlist::list.cbind() %>% base::as.data.frame()
      base::colnames(mean_exp) <- base::unique(hcobject[["data"]][[base::paste0("set", l, "_anno")]][[hcobject[["global_settings"]][["voi"]]]]) %>% base::as.character()
      l_co <- co[co %in% base::colnames(mean_exp)]
      final_colnames <<- c(final_colnames, l_co)
      
      mean_exp <- mean_exp[, l_co]
      
      return(mean_exp)
    }) %>% rlist::list.cbind()

    return(out)
  })
  title_plot <- base::names(hub_exp)[1]
  hub_exp <- hub_exp[[1]][stats::complete.cases(hub_exp[[1]]),]
  base::colnames(hub_exp) <- final_colnames
  
  
  p <- ComplexHeatmap::Heatmap(base::as.matrix(hub_exp), 
                               show_column_names = T, 
                               show_row_names = T, 
                               border = "black", 
                               col =  grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 7, name = "BrBG")))(50), 
                               cluster_columns = F, cluster_rows = F, rect_gp = grid::gpar(col = "black"),
                               heatmap_legend_param = list(title = ""), width = grid::unit(5, "cm"), height = grid::unit(7, "cm"),
                               row_names_gp = grid::gpar(fontsize = 8),
                               column_names_gp = grid::gpar(fontsize = 8), column_title = title_plot)
  return(p)
  
}

#' Change Colour Transparency
#' 
#' @param color String. Name of the colour.
#' @param alpha Integer between 1 and 100, giving the alpha value (saturation) of the target colour. Default is 100.
#' @noRd

makeTransparent <- function(color, alpha=100){

  newColor <- grDevices::col2rgb(color)
  transp <- base::apply(newColor, 2, function(rgbdata){grDevices::rgb(red=rgbdata[1], green=rgbdata[2],
                                              blue=rgbdata[3], alpha=alpha, maxColorValue=255)})
  return(transp)
}



#' Plot Boxplots From A List Of Dataframes As Input
#' @noRd

boxplot_from_list <- function(data, log_2, bool_plot){
  
  plts <- list()
  
  for(x in 1:base::length(data)){
    
    if(!base::is.data.frame(data[[x]])){
      
      print("Your list needs to contain data frames only!")
      
      break()
      
    }
    
    plts[[x]] <- boxplot_from_df(data[[x]], it = hcobject[["layers_names"]][x], log_2 = log_2, bool_plot = bool_plot)
    
  }
  return(plts)
}



#' Plot Boxplots From A Dataframe As Input
#' @noRd

boxplot_from_df <- function(data, it = NULL, log_2, bool_plot){
  
  # in case of large data set inform the user about the plots being split:
  if(base::ncol(data) > 50){
    
    print("To avoid over-crowded plots, several plots will be created with up to 50 samples each.")
    
  }
  
  if(log_2 == T){
    
    data[data < 1] <- NA
    
    data <- base::log(data, base = 2)
    
  }
  
  # if there are less than 50 samples, the data can be plotted in one plot:
  if(base::ncol(data) <= 50){
    
    bp_df <- base::data.frame(values= base::unlist(data), sample = base::rep(base::colnames(data), each = base::nrow(data)))
    
    p <- ggplot2::ggplot(bp_df, ggplot2::aes(x = sample, y = values)) + 
      ggplot2::geom_boxplot(fill = "gray90", color = "deepskyblue4")+
      ggplot2::theme_bw() + 
      ggplot2::ggtitle(it) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            text = ggplot2::element_text(size = 10))
    
    ggplot2::ggsave(filename = base::paste0("sample_distribution_bp", it,".pdf"), plot = p, device = cairo_pdf,
           path = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]],"/","sample_distribution_plots/"), width = 17, height = 7.8,  units = "in")
    
    return(p)
    
    # if there are more samples, the plots need to be split up to improve readability  
  }else{
    # find the best suitable number of boxplots per plot to make the figures as balanced as possible
    # (e.g. when there are 51 columns, we would not want a figure with 50 boxplots and another figure with just one boxplot
    # but rather one with 26 boxplots and one with 25)
    n <- find_best_mod(nc = base::ncol(data), ref = 50)
    
    plts <- list()
    # plot each plot: 
    for (x in 1: ((base::ncol(data) %/% n) + 1) ){
      
      start = ((x - 1) * n) + 1
      
      end = x * n
      
      if(start > base::ncol(data)){
        
        break()
        
      }
      
      if(end > base::ncol(data)){
        
        end <- base::ncol(data)
        
      }
      
      bp_df <- base::data.frame(values = base::unlist(data[, start:end]), sample = base::rep(base::colnames(data[, start:end]), each = base::nrow(data)))
      
      p <- ggplot2::ggplot(bp_df, ggplot2::aes(x = sample, y = values)) + 
        ggplot2::geom_boxplot(fill = "gray90", color = "deepskyblue4")+
        ggplot2::theme_bw() + 
        ggplot2::ggtitle(it) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
              text = ggplot2::element_text(size = 10))
      
      ggplot2::ggsave(filename = base::paste0("sample_distribution_bp", it, "_", x,".pdf"), plot = p, device = cairo_pdf,
             path = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]],"/","sample_distribution_plots/"), width = 10, height = 8,  units = "in")
      
      plts[[x]] <- p
      
    }
    return(plts)
  }
}



#' Plot Frequency Distributions From A List Of Dataframes As Input
#' @noRd

freqdist_plot_from_list <- function(data, log_2, bool_plot){
  
  plts <- list()
  
  for(x in 1:base::length(data)){
    
    if(!base::is.data.frame(data[[x]])){
      
      print("Your list needs to contain data frames only!")
      
      break()
      
    }
    
    plts[[x]] <- freqdist_plot_from_df(data[[x]], log_2 = log_2, it =  hcobject[["layers_names"]][x], bool_plot = bool_plot)
    
  }
  
  return(plts)
}




#' Plot Frequency Distributions From A Dataframe As Input
#' @noRd

freqdist_plot_from_df <- function(data, log_2, bool_plot, it = NULL){
  
  if(base::ncol(data) > 42){
    
    print("To avoid over-crowded plots, several plots will be created with up to 42 samples each.")
    
  }
  
  if(log_2 == T){
    
    data[data < 1] <- NA
    
    data <- base::log(data, base = 2)
    
  }
  
  if(base::ncol(data) <= 42){
    
    plot_list <- list()
    
    for(x in base::c(1:base::ncol(data))){
      
      plot_list[[base::paste0("plot_", x)]] <- ggplot2::ggplot(data, ggplot2::aes(x = data[, x])) + 
        ggplot2::geom_density(ggplot2::aes(y = ..count..), fill = "lightgray") +
        ggplot2::geom_vline(ggplot2::aes(xintercept = base::mean(data[,x][stats::complete.cases(data[, x])])), 
                   linetype = "dashed", size = 0.6,
                   color = "#FC4E07")+
        ggplot2::xlab(base::colnames(data)[x])+
        ggplot2::xlim(base::c(0, base::max(data[,x])))+
        ggplot2::theme_bw() + 
        ggplot2::theme(text = ggplot2::element_text(size = 10))
      
    }
    
    p <- base::do.call(gridExtra::grid.arrange, plot_list)
    
    ggplot2::ggsave(filename = base::paste0("sample_distribution_freq", it,".pdf"), plot = p, device = cairo_pdf,
           path = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]],"/","sample_distribution_plots/"), width = 10, height = 8, units = "in")
    
    return(p)
    
  }else{
    
    n <- find_best_mod(nc = base::ncol(data), ref = 42)
    
    plts <- list()
    
    for (x in 1: ((base::ncol(data) %/% n) + 1) ){
      
      start = ((x - 1) * n) + 1
      
      end = x * n
      
      if(start > base::ncol(data)){
        
        break()
        
      }
      
      if(end > base::ncol(data)){
        
        end <- base::ncol(data)
        
      }
      
      plot_list <- list()
      
      for(y in start:end){
        
        plot_list[[base::paste0("plot_",y)]] <- ggplot2::ggplot(data, ggplot2::aes(x = data[,y])) + 
          ggplot2::geom_density(ggplot2::aes(y = ..count..), fill = "lightgray") +
          ggplot2::geom_vline(ggplot2::aes(xintercept = base::mean(data[,y][stats::complete.cases(data[,y])])), 
                     linetype = "dashed", size = 0.6,
                     color = "#FC4E07")+
          ggplot2::xlab(base::colnames(data)[y])+
          ggplot2::xlim(base::c(0, base::max(data[,y])))+
          ggplot2::theme_bw() + 
          ggplot2::theme(text = ggplot2::element_text(size = 10))
        
      }
      
      p <- base::do.call(gridExtra::grid.arrange, plot_list)
      
      ggplot2::ggsave(filename = base::paste0("sample_distribution_freq", it, "_", x, ".pdf"), plot = p, device = cairo_pdf,
             path = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]],"/","sample_distribution_plots/"), width = 17, height = 7.8,  units = "in")
      
      plts[[x]] <- p
    }
    
    return(plts)
  }
}



#' Plot List Of Plots
#' @noRd

plot_list_of_plots <- function(plts){
  
  for (j in 1:base::length(plts)){
    
    tmp <- plts[[j]]
    
    if(ggplot2::is.ggplot(tmp)){
      
      graphics::plot(tmp)
      
    }else{
      
      for(p in tmp){
        
        graphics::plot(p)
        
      }
    }
  }
}



#' Split A Number In Well-Balanced Way
#' @noRd

find_best_mod <- function(nc, ref){
  
  mods <- NULL
  
  for(x in 1:ref){
    
    mods <- dplyr::bind_rows(mods, base::data.frame(x = x, mod = nc %% x))
    
  }
  
  mods <- mods[!base::duplicated(mods$mod),]
  
  mods <- mods[mods$mod == base::max(mods$mod),"x"][1]
  
  return(mods)
}



#' Filter Integrated Network For Non-White Genes
#' @noRd

network_filt <- function(){

  base::set.seed(123)
  gtc <- GeneToCluster()
  network <- hcobject[["integrated_output"]][["merged_net"]]
  del_v <- igraph::V(network)$name[!igraph::V(network)$name %in% gtc$gene]
  network <- igraph::delete.vertices(network, del_v)
  return(network)

}


#' For each sample calculates the mean expression per cluster
#' @param set An integer. Number of the dataset, for which samples the means should be calculated.
#' @noRd

sample_wise_cluster_expression <- function(set){

  gtc <- GeneToCluster()
  counts <- hcobject[["data"]][[base::paste0("set", set, "_counts")]]
  cluster_means <- base::lapply(base::unique(gtc$color[!gtc$color == "white"]), function(c){
    genes <- dplyr::filter(gtc, color == c) %>% dplyr::pull(., "gene")
    tmp <- counts[genes, ]
    tmp <- tmp[stats::complete.cases(tmp), ] %>% base::apply(., 2, base::mean)
  }) %>% rlist::list.rbind() %>% base::as.data.frame()
  base::rownames(cluster_means) <- base::unique(gtc$color[!gtc$color == "white"])

  return(cluster_means)

}


#' Hub Node Detection
#' 
#' Subroutine to find_hubs().
#' @noRd


hub_node_detection <- function(cluster, top = 10, save = T, tree_layout = F, TF_only = F, Plot){
  
  
  # extract chosen cluster as an isolated network:
  g <- cluster_to_network(cluster = cluster)
  # determine hub nodes:
  hub_out <- get_hub_nodes(network = g, top = top, TF_only = TF_only)
  
  # rank_df <- hub_out$rank_df
  if(base::nrow(hub_out$rank_df) == 0){
    return(hub_out)
  }

  # vertex size (10 for hub, 3 for non-hub):
  vertex_size <- base::lapply(igraph::V(g)$name, function(node){
    if(node %in% hub_out$hub_nodes){
      6
    }else{
      3
    } 
  }) %>% base::unlist()
  # labels (for hubs, none for non-hubs):
  label <- base::lapply(igraph::V(g)$name, function(node){
    if(node %in% hub_out$hub_nodes){
      node
    }else{
      NA
    } 
  }) %>% base::unlist()
  # layout:
  if(cluster == "all"){
    l <- hcobject[["integrated_output"]][["cluster_calc"]][["layout"]]
  }else{
    if(tree_layout){
      l <- igraph::layout_as_tree(g)
    }else{
      l <- igraph::layout.lgl(g) # issues with >1 graph components?
    }
    base::rownames(l) <- igraph::V(g)$name
  }

  
  # set some new node attributes:
  igraph::V(g)$size = vertex_size
  igraph::V(g)$label = NA
  igraph::V(g)$color = hub_out$colour_df$colour
  # plot network:
  network_with_labels(network = g, gene_labels = hub_out$hub_nodes, 
                      gene_ranks = base::c(1:base::length(hub_out$hub_nodes)),
                      l = l, label_offset = 10, 
                      title = base::c(cluster, top), save = save, Plot = Plot)
  

  
  return(hub_out)
}



#' Cluster To Network
#' 
#' Extracts a cluster as a standalone network.
#' @noRd

cluster_to_network <- function(cluster){

  gtc <- GeneToCluster() %>% dplyr::filter(., color == cluster)
  g <- hcobject[["integrated_output"]][["merged_net"]]
  g <- igraph::delete_vertices(g, igraph::V(g)$name[!igraph::V(g)$name %in% gtc$gene])
  
  return(g)

}




#' Get Hub Nodes
#' 
#' Subroutine to find_hubs().
#' @noRd

get_hub_nodes <- function(network = hcobject[["integrated_output"]][["merged_net"]], top = 10, TF_only){
  
  rank_df <- combined_centrality(network = network)
  rank_df$node <- base::rownames(rank_df)

  if(TF_only == "all"){
    cn <- base::grep(hcobject[["global_settings"]][["organism"]], base::colnames(hcobject[["supplementary_data"]][["TF"]]), ignore.case=TRUE, value=TRUE)
    rank_df <- dplyr::filter(rank_df, node %in% hcobject[["supplementary_data"]][["TF"]][[cn]])
  }
  if(!TF_only == "all" & !TF_only == FALSE){
    cn <- base::grep(hcobject[["global_settings"]][["organism"]], base::colnames(hcobject[["supplementary_data"]][["TF"]]), ignore.case=TRUE, value=TRUE)
    tmp <- hcobject[["supplementary_data"]][["TF"]][hcobject[["supplementary_data"]][["TF"]][, base::ncol(hcobject[["supplementary_data"]][["TF"]])] == TF_only,]
    rank_df <- dplyr::filter(rank_df, node %in% tmp[[cn]])
  }
  
  if(top > base::nrow(rank_df)){
    hub_nodes <- rank_df %>% dplyr::pull(., "node")
  }else{
    hub_nodes <- rank_df[1:top,] %>% dplyr::pull(., "node")
  }
  
  if(base::nrow(rank_df) == 0){
    message("No hubs found for this cluster.")
    return(list(rank_df = rank_df, colour_df = NULL, hub_nodes = hub_nodes))
  }
  colour_df <- centrality_colours(rank_df = rank_df, network = network)
  
  return(list(rank_df = rank_df, colour_df = colour_df, hub_nodes = hub_nodes))
}



#' Combined Centrality Measure
#' 
#' Subroutine to find_hubs().
#' @noRd

combined_centrality <- function(network){
  dc <- weighted_DC(network)
  cc <- weighted_CC(network)
  bc <- weighted_BC(network)
  
  rank_df <- base::data.frame(dc = base::rank(dc, ties.method = "average"),
                        cc = base::rank(cc, ties.method = "average"),
                        bc = base::rank(bc, ties.method = "average")
                        )
  
  rank_df$sum <-base::apply(rank_df, 1, base::sum)
  rank_df$id <- igraph::V(network) %>% base::as.character()
  
  # order based on highest rank (strongest hub candidates):
  rank_df <- rank_df[base::order(rank_df$sum, decreasing = T),]
  return(rank_df)
}


#' Weighted Degree Centrality 
#' 
#' Caclulates the weighted degree centrality of a node. Modified to be weighted from https://doi.org/10.1155/2019/9728742.
#' @noRd

weighted_DC <- function(network){
  message("Calculating weighted degree centrality.")
  out <- base::lapply(igraph::V(network), function(node){
            # get edges of node:
            e <- igraph::incident(graph = network, v = node, mode = "all")
            # get weights of edges:
            w <- igraph::edge_attr(network, "weight", index = e)
            # get all edge weights in the network:
            wall <- igraph::edge_attr(network, "weight", index = igraph::E(network))
            # return weighted degree centrality:
            return(base::sum(w)/base::sum(wall))
          }) %>% base::unlist()
  return(out)
}



#' Weighted Closeness Centrality 
#' 
#' Caclulates the weighted closeness centrality of a node. Modified to be weighted from https://doi.org/10.1155/2019/9728742.
#' @noRd

weighted_CC <- function(network){
  message("Calculating weighted closeness centrality.")
  # get |v|x|v| distance table with weighted shortest distances between all nodes:
  if(igraph::count_components(network) > 1){
    
  }else{
    dt <- igraph::distances(network, 
                            mode = "all", 
                            weights = 1-igraph::edge_attr(network, "weight", index = igraph::E(network))) 
  }
  
  # get sum of shortest distances for each node:
  d_sum <- base::apply(dt,2,sum)
  # return weighted closeness centrality for each node:
  return(igraph::vcount(network)/d_sum)
}



#' Weighted Betweenness Centrality 
#' 
#' Caclulates the weighted betweenness centrality of a node. Modified to be weighted from https://doi.org/10.1155/2019/9728742.
#' @noRd

weighted_BC <- function(network){
  message("Calculating weighted betweenness centrality.")
  # find shortest path for all pairs of nodes:
  out <- base::lapply(igraph::V(network), function(node){
    igraph::all_shortest_paths(graph = network, 
                                     from = node, 
                                     to = igraph::V(network)[!igraph::V(network) == node], 
                                     mode = "all",  
                                     weights = 1-igraph::edge_attr(network, "weight", index = igraph::E(network))) %>% 
      magrittr::extract2(1)
  }) 
  # combine list of lists into on long list:
  out <- base::do.call(base::c, out)
  # total number of shortest paths:
  total_sp <- base::length(out)
  # initialize counter for each node, set value to be zero:
  counts <- base::rep(0, igraph::vcount(network))
  # count number of shortest paths that pass each node:
  base::lapply(out, function(x){
    counts[x] <<- counts[x] + 1
  })
  # return weighted betweenness centrality per node:
  return(counts/total_sp)
}




#' Colours Based On Centrality
#' 
#' Subroutine to find_hubs().
#' @noRd

centrality_colours <- function(rank_df, network){
  
  mypalette <- grDevices::colorRampPalette(base::c('#FFF7EC', '#FC8D59', '#7F0000'))
  colour_df <- base::data.frame(sum = rank_df$sum, 
                          id = rank_df$id,
                          colour = (mypalette(base::max(rank_df$sum))[rank_df$sum]))
  colour_df <- colour_df[base::match(igraph::V(network), colour_df$id),]
  return(colour_df)

}



#' Plots A Network With Labels
#' 
#' Subroutine to find_hubs().
#' @noRd

network_with_labels <- function(network, gene_labels, gene_ranks, l, label_offset, title, save, Plot){
  
  
  plot_title <- base::paste0("cluster '", title[1], "' hub genes [top ",title[2], "]")
  
  new_genes_df <- base::data.frame(name = base::paste0("label_", base::c(1:base::length(gene_labels))), label = gene_labels %>% base::as.character())
  mean_coord_x <- stats::median(l[,1])
  mean_coord_y <- stats::median(l[,2])
  new_indeces <- base::match(new_genes_df$label, igraph::get.vertex.attribute(network)$name)
  new_genes_df$coords_x <- l[new_indeces,1]
  new_genes_df$coords_y <- l[new_indeces,2]
  left_up <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y >= mean_coord_y)%>%
    dplyr::pull(., label)
  left_down <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y < mean_coord_y)%>%
    dplyr::pull(., label)
  right_up <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y >= mean_coord_y)%>%
    dplyr::pull(., label)
  right_down <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y < mean_coord_y)%>%
    dplyr::pull(., label)
  
  
  
  new_position_l <- base::matrix(base::cbind(base::rep(base::ceiling(base::min(l[,1]))-label_offset, base::length(left_up)+ base::length(left_down)), 
                                 base::seq(from=base::ceiling(base::max(l[,2])), to = base::ceiling(base::min(l[,2])), 
                                     length.out = base::length(left_up)+ base::length(left_down))), ncol = 2)
  new_position_r <- base::matrix(base::cbind(base::rep(base::ceiling(base::max(l[,1]))+label_offset, base::length(right_up)+ base::length(right_down)), 
                                 base::seq(from=base::ceiling(base::max(l[,2])), to = base::ceiling(base::min(l[,2])), 
                                     length.out = base::length(right_up)+ base::length(right_down))), ncol = 2)
  new_position <- base::rbind(new_position_l, new_position_r)
  
  base::colnames(new_position) <- base::colnames(l)
  l2 <- base::rbind(l, new_position)%>%base::as.matrix()
  new_genes_df_l <- new_genes_df[new_genes_df$label %in% base::c(left_up, left_down),]
  new_genes_df_l <- new_genes_df_l[base::order(new_genes_df_l$coords_y, decreasing = T),]
  new_genes_df_r <- new_genes_df[new_genes_df$label %in% base::c(right_up, right_down),]
  new_genes_df_r <- new_genes_df_r[base::order(new_genes_df_r$coords_y, decreasing = T),]
  new_genes_df <- base::rbind(new_genes_df_l, new_genes_df_r)
  
  
  
  network2 <- igraph::add.vertices(network, nv = base::length(new_genes_df$name), 
                                   attr = list(name = new_genes_df$name))
  new_labels <- base::lapply(igraph::get.vertex.attribute(network2)$name, function(x){
    if(x %in% igraph::get.vertex.attribute(network)$name){
      NA
    }else{
      tmp_label <- dplyr::filter(new_genes_df, name == x)%>%
        dplyr::pull(., "label")
      tmp_rank <- gene_ranks[base::match(tmp_label, gene_labels)]
      base::paste0(tmp_label, " [",tmp_rank, "]")
    }
  })%>% base::unlist()
  
  igraph::V(network2)$label <- new_labels

  
  name_to_id <- base::data.frame(name = igraph::V(network2)$name, id = 1:base::length(igraph::V(network2)$name))
  
  
  nti1 <- name_to_id[name_to_id$name %in% new_genes_df[,1],]
  base::rownames(nti1) <- nti1$name
  nti1 <- nti1[new_genes_df[,1],]
  
  nti2 <- name_to_id[name_to_id$name %in% new_genes_df[,2],]
  base::rownames(nti2) <- nti2$name
  nti2 <- nti2[new_genes_df[,2],]
  
  new_edges <- base::matrix(base::c(nti1$id, nti2$id), ncol = 2, byrow = F)
  
  new_edges <- base::as.vector(base::t(new_edges))

  network2 <- igraph::add.edges(network2, new_edges)
  
  new_edge_color <- base::apply(igraph::get.edgelist(network2),1, function(x){
    
    if(x[2] %in% base::as.character(new_genes_df$name)){
      
      "#2B8CBE"
    }else{
      "lightgrey"
    }
  })
  
  vertex_shape <- base::lapply(igraph::V(network2)$name, function(x){
    if(x %in% igraph::V(network)$name){
      "circle"
    }else{
      "none"
    }
  }) %>% base::unlist()
  
  new_label_color <- base::lapply(igraph::get.vertex.attribute(network2)$name, function(x){
    if(x %in% new_genes_df$name){
      tmp_gene_n <- dplyr::filter(new_genes_df, name == x)%>%
        dplyr::pull(., label)
      tfs <- hcobject[["supplementary_data"]][["TF"]][, base::grep(base::colnames(hcobject[["supplementary_data"]][["TF"]]), pattern = hcobject[["global_settings"]][["organism"]], ignore.case = T)]
      if(tmp_gene_n %in% tfs){
        "goldenrod"
      }else{
        igraph::get.vertex.attribute(network2, name = "color", index = tmp_gene_n)
      }
    }else{
      NA
    }
  })%>% base::unlist()
  
  name_to_size <- base::data.frame(name = igraph::V(network)$name, size = igraph::V(network)$size)
  
  vertex_size <- base::lapply(igraph::V(network2)$name, function(x){
    if(x %in% igraph::V(network)$name){
      dplyr::filter(name_to_size, name == x) %>% dplyr::pull(., "size")
    }
    else{
      0
    }
  }) %>% base::unlist()
  if(save == T){
    Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", base::paste0("cluster ", title[1], " hub genes"), ".pdf"), 
                    width = 20, height = 15)
    
    igraph::plot.igraph(network2, vertex.size = vertex_size, vertex.label = new_labels, vertex.label.cex = 1.5,
                        layout = l2, vertex.label.dist = 1, vertex.shape = vertex_shape,
                        edge.color = new_edge_color, vertex.label.color = new_label_color)
    graphics::title(plot_title, cex.main=3)
    
    grDevices::dev.off()
  }
  if(Plot){
    igraph::plot.igraph(network2, vertex.size = vertex_size, vertex.label = new_labels, vertex.label.cex = 0.75,
                        layout = l2, main = plot_title, vertex.label.dist = 1, vertex.shape = vertex_shape,
                        edge.color = new_edge_color, vertex.label.color = new_label_color)
  }
  
  
}




#' Plots Th Expression Of Hub Genes As Heatmaps
#' 
#' Subroutine to find_hubs().
#' @noRd

plot_hub_exp <- function(hubs_df = hcobject[["integrated_output"]][["hub_out"]]){
  
  # create save folder:
  if(!"hub_expression_plots" %in% base::list.dirs(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]],"/"))){
    
    base::dir.create(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]],"/","hub_expression_plots"))
    
  }
  
  
  #get column order from module heatmap:
  co <- ComplexHeatmap::column_order(hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]]) %>% base::unlist() %>% base::as.numeric()
  co <- base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]][,-base::ncol(hcobject[["integrated_output"]][["GFC_all_layers"]])])[co]
  hub_exp <- base::apply(hubs_df, 2, function(x){

    genes <- x[!x == " "]
 
    if(base::length(genes) != 0){
      out <- base::lapply(1:base::length(hcobject[["layers"]]), function(l){
        exp <- hcobject[["data"]][[base::paste0("set", l, "_counts")]][genes,]
        mean_exp <- base::lapply(base::unique(hcobject[["data"]][[base::paste0("set", l, "_anno")]][[hcobject[["global_settings"]][["voi"]]]]), function(c){
          samples <- base::rownames( hcobject[["data"]][[base::paste0("set", l, "_anno")]][hcobject[["data"]][[base::paste0("set", l, "_anno")]][[hcobject[["global_settings"]][["voi"]]]]==c,] )
          tmp <- dplyr::select(exp, samples) %>% base::rowMeans()
        }) %>% rlist::list.cbind() %>% base::as.data.frame()
        
        base::colnames(mean_exp) <- base::unique(hcobject[["data"]][[base::paste0("set", l, "_anno")]][[hcobject[["global_settings"]][["voi"]]]]) %>% base::as.character()

        # l_co <- co[co %in% base::colnames(mean_exp)]
        # mean_exp <- dplyr::select(mean_exp, tidyselect::all_of(l_co))

        return(mean_exp)
      }) #%>% rlist::list.cbind()
      
      return(out)
    }else{
      return(NULL)
    }
    
    
  })
  pl <- base::lapply(1:base::length(hub_exp), function(x){
    if(!base::is.null(hub_exp[[x]])){
      pl2 <- base::lapply(1:base::length(hub_exp[[x]]), function(y){
        
        if(!base::is.null(hub_exp[[x]][[y]])){
          p <- pheatmap::pheatmap(hub_exp[[x]][[y]], 
                                  cluster_cols = F, 
                                  cluster_rows = F, 
                                  border_color = "black", 
                                  cellwidth = 20, 
                                  cellheight = 20, 
                                  main = base::paste0(base::colnames(hubs_df)[x], " dataset ", y),
                                  fontsize_col = 14, 
                                  color = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 7, name = "BrBG")))(50),
                                  scale = "row")
          return(p$gtable)
        }
      })
      cp1 <- cowplot::plot_grid(plotlist = pl2, nrow = 1, align = "h")
      
      grDevices::png(filename = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]],"/hub_expression_plots", 
                                             "/Hub_TF_expr_", base::colnames(hubs_df)[x], ".png"), 
                     res = 300, width = 3500, height = 2500)
      
      graphics::plot(cp1)
      
      grDevices::dev.off()
      
    }
    
  })
  

  
}


#' Run All Clustering Algorithms
#' 
#' Runs all clustering algorithms on the network and returns the clustering.
#' @noRd

run_all_cluster_algos <- function(){

    output <- list()

    network <- hcobject[["integrated_output"]][["merged_net"]]

    cluster_algo_list <- base::c("cluster_louvain",
                         "cluster_fast_greedy",
                         "cluster_infomap",
                         "cluster_walktrap",
                         "cluster_label_prop",
                         "cluster_leiden")


    color.cluster <- get_cluster_colours()


    current_algo <- NULL

    for(c in base::unique(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]$color)){
      genes <- dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], color == c) %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)
      current_algo <- base::rbind(current_algo, base::data.frame(gene = genes, cluster = base::rep(c, base::length(genes))))
    }

    current_algo$cluster <- base::as.factor(current_algo$cluster)
    
    output[[hcobject[["global_settings"]][["chosen_clustering_algo"]]]] <- current_algo

    for(a in cluster_algo_list[!cluster_algo_list == hcobject[["global_settings"]][["chosen_clustering_algo"]]]){

    base::set.seed(1)
    base::set.seed(.Random.seed[1])

    # add exception for leiden:
    if(a == "cluster_leiden"){
      partition <- leidenAlg::leiden.community(graph = network)
      partition_df <- base::data.frame(gene = partition$names, cluster = base::as.numeric(partition$membership))
      partition_df$cluster <- partition_df$cluster + 1
    }else{
      cfg = base::get(a)(network)
      partition_df <- base::data.frame(gene = igraph::get.vertex.attribute(network)$name, cluster = cfg$membership)
      
    }

    partition_df$cluster <- base::lapply(partition_df$cluster, function(x){
      color.cluster[x]
    })%>% base::unlist()
    
    
    
    partition_df_table <- base::table(partition_df$cluster) %>% 
      base::as.data.frame()

    partition_df <- partition_df[partition_df$cluster %in% (dplyr::filter(partition_df_table, Freq >= hcobject[["global_settings"]][["min_nodes_number_for_cluster"]]) %>%
                 dplyr::pull(., "Var1")),]

    partition_df$cluster <- base::as.factor(partition_df$cluster)

    output[[a]] <- partition_df
    }
    return(output)
}


#' Plot PCA of top most variant
#' 
#' Subroutine to PCA_algo_compare().
#' @noRd

plot_PCA_topvar <- function(PCA_save_folder){

  plotlist <- list()
  pca_list <- list()
  for(i in 1:base::length(hcobject[["layers"]])){

    pca <- stats::prcomp(base::t(hcobject[["layer_specific_outputs"]][[base::paste0("set",i)]][["part1"]][["topvar"]]), scale = T)
    pca.var <- pca$sdev^2
    pca.var.per <- base::data.frame(pc = 1:base::length(pca.var), val = base::round(pca.var/base::sum(pca.var)*100,1))
    pca.data <- base::data.frame(Sample = base::rownames(pca$x),
                           X = pca$x[,1],
                           Y = pca$x[,2],
                           Group = hcobject[["data"]][[base::paste0("set",i, "_anno")]][[hcobject[["global_settings"]][["voi"]]]])

    num_colours <- base::length(base::unique(dplyr::pull(hcobject[["data"]][[base::paste0("set",i,"_anno")]], hcobject[["global_settings"]][["voi"]])))
    my_palette <- ggsci::pal_d3("category20")(20)
    if(base::length(num_colours) > base::length(my_palette)){
      my_colours <- grDevices::colorRampPalette(my_palette)(num_colours)
    }else{
      my_colours <- my_palette[1:num_colours]
    }
    

    p <- ggplot2::ggplot(pca.data, ggplot2::aes(x = X, y = Y, col = Group, label = Sample))+
      ggplot2::geom_point(size = 4)+
      ggplot2::ylab(base::paste0("PC 2", " (", pca.var.per[2, 2], "%)"))+
      ggplot2::xlab(base::paste0("PC 1", " (", pca.var.per[1, 2], "%)"))+
      ggplot2::theme_bw() + 
      ggplot2::ggtitle(base::paste0(hcobject[["layers_names"]][i], " by topvar"))+
      ggplot2::scale_color_manual(values = my_colours)


      folder <- PCA_save_folder
      if(!folder %in% base::list.dirs(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]]))) {

        base::dir.create(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/",  folder))

      }

    Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", folder, "/PCA_topvar_", hcobject[["layers_names"]][i], ".pdf"), width = 10, height = 7)
    graphics::plot(p)
    grDevices::dev.off()
    plotlist[[i]] <- p
    pca_list[[i]] <- pca
  }
  if(base::length(hcobject[["layers"]]) == 1){
    cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 1)
  }else{
    cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 2, align = "hv")
  }
  graphics::plot(cp)
  return(pca_list)

}


#' Plot PCA based on cluster expressions
#' 
#' Subroutine to PCA_algo_compare().
#' @noRd

plot_PCA_cluster <- function(gtc = NULL, algo = NULL, PCA_save_folder){

  plotlist <- list()
  pca_list <- list()
  for(l in 1:base::length(hcobject[["layers"]])){
    
    if(base::is.null(gtc)){
      FC <- intra_sample_FC(l)
    }else{
      FC <- intra_sample_FC(l, gtc = gtc)

    }

    pca <- stats::prcomp(base::t(FC), scale = T)
    pca.var <- pca$sdev^2
    pca.var.per <- base::data.frame(pc = 1:base::length(pca.var), val = base::round(pca.var/base::sum(pca.var)*100,1))
    anno <- hcobject[["data"]][[base::paste0("set", l, "_anno")]]
    if(base::ncol(pca$x) == 1){
      pca.data <- base::data.frame(Sample = base::rownames(pca$x),
                             X = pca$x[,1],
                             Y = 0,
                             Group = dplyr::pull(anno , hcobject[["global_settings"]][["voi"]]))
    }else{
      pca.data <- base::data.frame(Sample = base::rownames(pca$x),
                             X = pca$x[,1],
                             Y = pca$x[,2],
                             Group = dplyr::pull(anno , hcobject[["global_settings"]][["voi"]]))
    }
    num_colours <- base::length(base::unique(dplyr::pull(hcobject[["data"]][[base::paste0("set",l,"_anno")]], hcobject[["global_settings"]][["voi"]])))
    my_palette <- ggsci::pal_d3("category20")(20)
    if(base::length(num_colours) > base::length(my_palette)){
      my_colours <- grDevices::colorRampPalette(my_palette)(num_colours)
    }else{
      my_colours <- my_palette[1:num_colours]
    }

    p <- ggplot2::ggplot(pca.data, ggplot2::aes(x = X, y = Y,  label= Sample, color = Group))+
      ggplot2::geom_point(size = 4)+
      ggplot2::scale_color_manual(values = my_colours)+
      ggplot2::theme_bw() + 
      ggplot2::ggtitle(base::paste0(hcobject[["layers_names"]][l], " by module - ", algo))+
      ggplot2::ylab(base::paste0("PC 2", " (", pca.var.per[2, 2], "%)"))+
      ggplot2::xlab(base::paste0("PC 1", " (", pca.var.per[1, 2], "%)"))
    plotlist[[l]] <- p
    pca_list[[l]] <- pca

    folder <- PCA_save_folder
    if(!folder %in% base::list.dirs(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]]))) {

      base::dir.create(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", folder))

    }

    Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/",  folder, "/PCA_module_", hcobject[["layers_names"]][l], "_", algo, ".pdf"), width = 10, height = 7)
    graphics::plot(p)
    grDevices::dev.off()
  }
  if(base::length(hcobject[["layers"]]) == 1){
    cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 1)
  }else{
    cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 2, align = "hv")
  }
 graphics::plot(cp)
 return(pca_list)
}



#' Calculate Intra-Sample Fold Changes
#' 
#' Mean cluster exression from mean sample expression.
#' Subroutine to PCA_algo_compare().
#' @noRd

intra_sample_FC <- function(l, gtc = NULL){

  if(base::is.null(gtc)){
    gene_to_cluster <- GeneToCluster()
  }else{
    base::colnames(gtc) <- base::c("gene", "cluster")
    gene_to_cluster <- gtc
    gene_to_cluster$gene <- base::as.character(gene_to_cluster$gene)
    gene_to_cluster$cluster <- base::as.character(gene_to_cluster$cluster)
    base::colnames(gene_to_cluster) <- base::c("gene", "color")
  }


  counts <- hcobject[["data"]][[base::paste0("set", l , "_counts")]]

  mean_expression_per_sample <- base::apply(counts, 2, base::mean)

  mean_expression_per_cluster <- NULL

  for(c in base::unique(gene_to_cluster$color)){
    if(!c == "white"){
      genes <- gene_to_cluster[gene_to_cluster$color == c, ] %>%
        dplyr::pull(., "gene")

      filt_counts <- counts[base::rownames(counts) %in% genes, ]
      tmp <- base::apply(filt_counts, 2, base::mean)%>%
        base::as.data.frame()%>%
        base::t()%>%
        base::as.data.frame()
      base::rownames(tmp) <- c
      base::colnames(tmp) <- base::colnames(counts)
      mean_expression_per_cluster <- base::rbind(mean_expression_per_cluster, tmp)
    }
  }
  mean_expression_per_cluster <- base::rbind(mean_expression_per_cluster, mean_expression_per_sample)
  FC_from_mean_per_cluster <- base::apply(mean_expression_per_cluster,2, function(x){
    x/x[base::length(x)]
  }) %>% base::as.data.frame()
  base::colnames(FC_from_mean_per_cluster) <- base::colnames(mean_expression_per_cluster)
  base::rownames(FC_from_mean_per_cluster) <- base::rownames(mean_expression_per_cluster)
  FC_from_mean_per_cluster <- FC_from_mean_per_cluster[1:(base::nrow(FC_from_mean_per_cluster)-1),]

  return(FC_from_mean_per_cluster)

}



#' Find New Control
#' 
#' Detects sample subset with highest control content. Subroutine to cut_hclust().
#' @noRd

find_new_ctrl <- function(anno, l){
  ctrl_samples <- anno[base::grepl(hcobject[["global_settings"]][["control"]], anno[[base::paste0(hcobject[["global_settings"]][["voi"]], "_old")]], ignore.case = T),] %>% base::rownames()
  new_cons <- base::unique(dplyr::pull(anno, hcobject[["global_settings"]][["voi"]]))
  
  maxvec <- base::lapply(new_cons, function(x){
    tmp <- anno[anno[[hcobject[["global_settings"]][["voi"]]]] == x,]
    l <- base::intersect(ctrl_samples, base::rownames(tmp)) %>% base::length()
    return(l)
  }) %>% base::unlist()
  new_ctrl <- new_cons[base::which(maxvec == base::max(maxvec))[1]]
  base::print(base::paste0("New control: ", new_ctrl))
  anno[hcobject[["global_settings"]][["voi"]]][anno[hcobject[["global_settings"]][["voi"]]] == new_ctrl] <- base::paste0(hcobject[["layers_names"]][l], "_",  hcobject[["global_settings"]][["control"]])
  return(anno)
}


#' Get Annotation Matrix
#' 
#' Subroutine to col_anno_categorical().
#' @noRd

get_anno_matrix <- function(variables){
  all_mats <- list()
  for(v in 1:base::length(variables)){
    if(base::is.na(variables[v])){
      names <- dplyr::pull(hcobject[["data"]][[base::paste0("set", v, "_anno")]], hcobject[["global_settings"]][["voi"]]) %>% base::unique()
      layer_mat <- base::matrix(base::rep(0,base::length(names[names %in% base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])])), ncol = 1)
      base::rownames(layer_mat) <- names[names %in% base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])]
    }else{
      anno <- hcobject[["data"]][[base::paste0("set", v, "_anno")]]
      if(!variables[v] %in% colnames(anno)){
        stop(variables[v], " is not a column name found in the annotation of dataset ", v, ". Please check the spelling.")
      }
      layer_mat <- NULL
      for(i in base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])[!base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]]) == "Gene"]){

        if(i %in% dplyr::pull(anno, hcobject[["global_settings"]][["voi"]])){
          tmp <- dplyr::filter(anno, anno[hcobject[["global_settings"]][["voi"]]] == i)%>%
            dplyr::pull(., variables[v])

          tmp <- base::table(tmp)%>%
            base::t()%>%
            base::as.data.frame()

          tmp$Var1 <- NULL
          base::colnames(tmp) <- c("var", "freq")
          mat <- base::matrix(tmp$freq, nrow = 1, byrow = T)
          base::colnames(mat) <- tmp$var
          base::rownames(mat) <- i
          layer_mat <- dplyr::bind_rows(base::as.data.frame(layer_mat),  base::as.data.frame(mat))
        }
      }
    }
    layer_mat[base::is.na(layer_mat)] <- 0
    layer_mat <- base::as.matrix(layer_mat)
    all_mats[[v]] <- layer_mat
  }
  return(all_mats)
}



#' Unify Matrices
#' 
#' Subroutine to col_anno_categorical().
#' @noRd

unify_mats <- function(mat_list){
  merged_mat <- NULL
  new_rn <- NULL
  for(x in mat_list){
    if(!base::is.null(merged_mat)){
      if(base::all(base::colnames(merged_mat) %in% base::colnames(x))){
        new_rn <- base::c(new_rn, base::rownames(x))
        x <- x[, base::colnames(merged_mat)]
        merged_mat <- base::rbind(merged_mat,x)
      }else{
        new_rn <- base::c(new_rn, base::rownames(x))
        merged_mat <- base::merge(merged_mat, x, by = "row.names", all = T)
        merged_mat$Row.names <- NULL
      }
    }else{
      new_rn <- base::c(new_rn, base::rownames(x))
      merged_mat <- base::merge(merged_mat, x, by = "row.names", all = T)
      merged_mat$Row.names <- NULL
    }

  }
  merged_mat[base::is.na(merged_mat)] <- 0
  merged_mat <- merged_mat[, base::colSums(merged_mat)>0]
  base::rownames(merged_mat) <- new_rn
  return(merged_mat)
}




#' Adapted Version Of plot_cluster_heatmap()
#' 
#' Subroutine to regroup_cluster_heatmap().
#' @noRd

replot_cluster_heatmap <- function(col_order = NULL, 
                                       row_order = NULL, 
                                       cluster_columns = T,
                                       cluster_rows = T, 
                                       k = 0, 
                                       return_HM = T, 
                                       cat_as_bp = NULL, 
                                       file_name = "module_heatmap.pdf",
                                       GFCs,
                                       group,
                                       data){


  # set user specific enrichments if they exist:
  if("enriched_per_cluster" %in% base::names(hcobject[["satellite_outputs"]])){
    if(!base::is.null(hcobject[["satellite_outputs"]][["enriched_per_cluster"]])){
      user_enrichment_1 <- hcobject[["satellite_outputs"]][["enriched_per_cluster"]][["categories_per_cluster"]]
    }
  }else{
    user_enrichment_1 <- NULL
  }
  
  if("enriched_per_cluster2" %in% base::names(hcobject[["satellite_outputs"]])){
    if(!base::is.null(hcobject[["satellite_outputs"]][["enriched_per_cluster2"]])){
      user_enrichment_2 <- hcobject[["satellite_outputs"]][["enriched_per_cluster2"]][["categories_per_cluster"]]
    }
  }else{
    user_enrichment_2 <- NULL
  }

  column_anno_categorical <- NULL

  column_anno_numerical <- NULL
  
  if(base::is.null(cat_as_bp)){
    if(!base::is.null(column_anno_categorical)){
      cat_as_bp <- base::rep(F, base::length(column_anno_categorical))
    }
  }
  
  base::gc()

  # filter for included clusters (non-white)
  c_df <- dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], cluster_included == "yes")
  mat_heatmap <- NULL


  if(!base::is.null(row_order)){
    for (c in row_order){
      #get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)

      # GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(GFCs, Gene %in% genes)
      c_GFC_means <- base::apply(c_GFCs[, base::c(1:(base::ncol(c_GFCs)-1))], 2, base::mean)

      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)

    }
    base::rownames(mat_heatmap) <- row_order

  }else{
    for (c in base::unique(c_df$color)){
      #get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)


      # GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(GFCs, Gene %in% genes)

      if(base::is.vector(c_GFCs)){
        c_GFC_means <- cGFCs
      }else{
        c_GFC_means <- base::apply(c_GFCs[, base::c(1:(base::ncol(c_GFCs)-1))] %>% 
          base::as.data.frame(), 2, base::mean)
      }


      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)

    }
    base::rownames(mat_heatmap) <- c_df$color
  }


  base::colnames(mat_heatmap) <- base::colnames(GFCs)[1:(base::ncol(GFCs)-1)]

  if(!base::is.null(col_order)){
    mat_heatmap <- mat_heatmap %>% base::as.data.frame()
    mat_heatmap <- mat_heatmap[, col_order] %>% base::as.matrix()
  }

  enrich_mat1 <- list()
  enrich_count1 <- list()
  enrich_mat2 <- list()
  enrich_count2 <- list()

  if(!base::is.null(row_order)){
    if(!base::is.null(user_enrichment_1)){
      for(x in row_order){
        enrich_mat1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
          dplyr::pull(., count)
        enrich_count1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
          dplyr::pull(., hits) %>%
          dplyr::first(.)
      }
    }
    if(!base::is.null(user_enrichment_2)){
      for(x in row_order){
        enrich_mat2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
          dplyr::pull(., count)
        enrich_count2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
          dplyr::pull(., hits) %>%
          dplyr::first(.)
      }
    }

  }else{
    for(x in base::unique(user_enrichment_1$cluster)){
      enrich_mat1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
        dplyr::pull(., count)
      enrich_count1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
        dplyr::pull(., hits) %>%
        dplyr::first(.)
    }
    for(x in base::unique(user_enrichment_2$cluster)){
      enrich_mat2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
        dplyr::pull(., count)
      enrich_count2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
        dplyr::pull(., hits) %>%
        dplyr::first(.)
    }
  }

  if(!base::length(enrich_mat1) == 0){
    enrich_mat1 <- base::matrix(base::unlist(enrich_mat1), nrow = base::length(enrich_mat1), byrow = T)
    enrich_count1 <- base::unlist(enrich_count1)

  }
  if(base::all(enrich_mat1 == 0) == T){
    enrich_mat1 <- list()
  }

  if(!base::length(enrich_mat2) == 0){

    enrich_mat2 <- base::matrix(base::unlist(enrich_mat2), nrow = base::length(enrich_mat2), byrow = T)
    enrich_count2 <- base::unlist(enrich_count2)
  }
  if(base::all(enrich_mat2 == 0) == T){
    enrich_mat2 <- list()
  }

  if(!base::is.null(row_order)){
    cluster_colors <- base::factor(row_order)
    base::names(cluster_colors) <- row_order
    c_df <- c_df[base::match(row_order, c_df$color),]
  }else{
    cluster_colors <- base::factor(c_df$color)
    base::names(cluster_colors) <- c_df$color
    row_order <- base::unique(c_df$color)
  }


  if(base::length(enrich_mat1) == 0 & base::length(enrich_mat2) == 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = ComplexHeatmap::anno_simple(row_order, col = cluster_colors,
                                                                   simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")),
                                            genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(2.5, "cm")),
                                            gene_nums = ComplexHeatmap::anno_text(c_df$gene_no, width = grid::unit(1.5, "cm"), gp = grid::gpar(fontsize = 10)),


                                            which = "row",
                                            width = grid::unit(4.5, "cm"),
                                            annotation_name_side = "top",
                                            gap = grid::unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = grid::gpar(fontsize = 8))

    lgd_list <- list(

    )
  }
  else if(base::length(enrich_mat1) > 0 & base::length(enrich_mat2) == 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = ComplexHeatmap::anno_simple(row_order, col = cluster_colors,
                                                                   simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")),
                                            genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(2.5, "cm")),

                                            enriched_count = ComplexHeatmap::anno_text(base::paste0(enrich_count1, "/", c_df$gene_no), width = grid::unit(1.5, "cm")),
                                            enriched = ComplexHeatmap::anno_barplot(enrich_mat1,
                                                                      width = grid::unit(3, "cm"),
                                                                      gp = grid::gpar(fill = RColorBrewer::brewer.pal(n = 12, name = "Paired"),
                                                                                col = RColorBrewer::brewer.pal(n = 12, name = "Paired"))),
                                            which = "row",
                                            width = grid::unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = grid::unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = grid::gpar(fontsize = 8))

    lgd_list <- list(

      ComplexHeatmap::Legend(labels = base::unique(user_enrichment_1$cell_type), title = "enriched",
                             legend_gp = grid::gpar(col = RColorBrewer::brewer.pal(n = 12, name = "Paired")),
                             type = "points", pch = 15)
    )
  }
  else if(base::length(enrich_mat1) == 0 & base::length(enrich_mat2) > 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = ComplexHeatmap::anno_simple(row_order, col = cluster_colors,
                                                                   simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")),
                                            genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(1.5, "cm")),

                                            enriched_count = ComplexHeatmap::anno_text(base::paste0(enrich_count2, "/", c_df$gene_no), width = grid::unit(1.5, "cm"),
                                                                       gp = grid::gpar(fontsize = 8)),
                                            enriched = ComplexHeatmap::anno_barplot(enrich_mat2,
                                                                      width = grid::unit(5, "cm"),
                                                                      gp = grid::gpar(fill = ggsci::pal_d3(palette = "category20")(base::ncol(enrich_mat2)),
                                                                                col = ggsci::pal_d3(palette = "category20")(base::ncol(enrich_mat2)))),

                                            which = "row",
                                            width = grid::unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = grid::unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = grid::gpar(fontsize = 8))

    lgd_list <- list(


      ComplexHeatmap::Legend(labels = base::unique(user_enrichment_2$cell_type), title = "enriched",
                             legend_gp = grid::gpar(col = ggsci::pal_d3(palette = "category20")(20)),
                             type = "points", pch = 15)
    )
  }else{
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = ComplexHeatmap::anno_simple(row_order, col = cluster_colors,
                                                                   simple_anno_size = grid::unit(0.25, "cm"), gp = grid::gpar(col = "black")),
                                            genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(0.75, "cm")),

                                            enriched_count_1 = ComplexHeatmap::anno_text(base::paste0(enrich_count1, "/", c_df$gene_no),
                                                                         width = grid::unit(0.75, "cm"),
                                                                         gp = grid::gpar(fontsize = 8)),
                                            enriched_1 = ComplexHeatmap::anno_barplot(enrich_mat1,
                                                                      width = grid::unit(2, "cm"),
                                                                      gp = grid::gpar(fill = RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:base::ncol(enrich_mat1)],
                                                                                col = RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:base::ncol(enrich_mat1)]),
                                                                      baseline = 0),
                                            enriched_count_2 = ComplexHeatmap::anno_text(base::paste0(enrich_count2, "/", c_df$gene_no),
                                                                         width = grid::unit(0.75, "cm"),
                                                                         gp = grid::gpar(fontsize = 8)),
                                            enriched_2 = ComplexHeatmap::anno_barplot(enrich_mat2,
                                                                      width = grid::unit(2, "cm"),
                                                                      gp = grid::gpar(fill = ggsci::pal_d3(palette = "category20")(20)[1:base::ncol(enrich_mat2)],
                                                                                col = ggsci::pal_d3(palette = "category20")(20)[1:base::ncol(enrich_mat2)]),
                                                                      baseline = 0),
                                            which = "row",
                                            width = grid::unit(12, "cm"),
                                            annotation_name_side = "top",
                                            gap = grid::unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = grid::gpar(fontsize = 8))

    lgd_list <- list(


      ComplexHeatmap::Legend(labels = base::unique(user_enrichment_1$cell_type), title = "enriched_1",
                             legend_gp = grid::gpar(col = ggsci::pal_d3(palette = "category20")(20)[1:base::ncol(enrich_mat1)]),
                             type = "points", pch = 15),
      ComplexHeatmap::Legend(labels = base::unique(user_enrichment_2$cell_type), title = "enriched_2",
                             legend_gp = grid::gpar(col = ggsci::pal_d3(palette = "category20")(20)[1:base::ncol(enrich_mat2)]),
                             type = "points", pch = 15)
    )
  }


  anno_list <- NULL


  if(!base::length(column_anno_categorical) == 0){
    for(a in 1:base::length(column_anno_categorical)){
      base::set.seed(a)
      tmp_colour <- grDevices::colorRampPalette(ggsci::pal_d3("category20")(20))(base::ncol(column_anno_categorical[[a]]))
      if(cat_as_bp[a] == T){
        column_anno_categorical[[a]][base::is.na(column_anno_categorical[[a]])] <- 0
        if(base::is.null(anno_list)){
          anno_list <- ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_barplot(column_anno_categorical[[a]]%>%base::as.matrix(),
                                                                                             width = grid::unit(2, "cm"),
                                                                                             gp = grid::gpar(fill = tmp_colour,
                                                                                                       col = tmp_colour)),
                                                                     which = "column",
                                                                     height = grid::unit(1, "cm"),
                                                                     annotation_name_side = "right",
                                                                     gap = grid::unit(2, "mm"),
                                                                     annotation_name_rot = 0,
                                                                     annotation_name_gp = grid::gpar(fontsize = 8),
                                                                     annotation_label = base::names(column_anno_categorical)[a])
        }else{
          anno_list <-  ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_barplot(column_anno_categorical[[a]]%>%base::as.matrix(),
                                                                                             width = grid::unit(2, "cm"),
                                                                                             gp = grid::gpar(fill = tmp_colour,
                                                                                                       col = tmp_colour)),
                                                                     which = "column",
                                                                     height = grid::unit(1, "cm"),
                                                                     annotation_name_side = "right",
                                                                     gap = grid::unit(2, "mm"),
                                                                     annotation_name_rot = 0,
                                                                     annotation_name_gp = grid::gpar(fontsize = 8),
                                                                     annotation_label = base::names(column_anno_categorical)[a]), direction = "vertical")
        }
        
        
        
      }else{
        if(base::is.null(anno_list)){
          anno_list <- ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_lines(column_anno_categorical[[a]] %>% base::as.matrix(),width = grid::unit(2, "cm"),
                                                                                           gp = grid::gpar(col = tmp_colour),
                                                                                           add_points = TRUE,
                                                                                           pt_gp = grid::gpar(col = tmp_colour), pch = 16),
                                                                     which = "column",
                                                                     height = grid::unit(1, "cm"),
                                                                     annotation_name_side = "right",
                                                                     gap = grid::unit(2, "mm"),
                                                                     annotation_name_rot = 0,
                                                                     annotation_name_gp = grid::gpar(fontsize = 8),
                                                                     annotation_label = base::names(column_anno_categorical)[a])
        }else{
          anno_list <-  ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_lines(column_anno_categorical[[a]] %>% base::as.matrix(),width = grid::unit(2, "cm"),
                                                                                           gp = grid::gpar(col = tmp_colour),
                                                                                           add_points = TRUE,
                                                                                           pt_gp = grid::gpar(col = tmp_colour), pch = 16),
                                                                     which = "column",
                                                                     height = grid::unit(1, "cm"),
                                                                     annotation_name_side = "right",
                                                                     gap = grid::unit(2, "mm"),
                                                                     annotation_name_rot = 0,
                                                                     annotation_name_gp = grid::gpar(fontsize = 8),
                                                                     annotation_label = base::names(column_anno_categorical)[a]), direction = "vertical")
        }
        
      }
      
      

      lgd_list <- rlist::list.append(lgd_list, ComplexHeatmap::Legend(labels = colnames(column_anno_categorical[[a]]%>%as.matrix()), title = names(column_anno_categorical)[a],
                                                                      legend_gp = gpar(col = tmp_colour),
                                                                      type = "points", pch = 15))
    }

  }



  if(!base::length(column_anno_numerical) == 0){
    for(a in 1:base::length(column_anno_numerical)){
      tmp_col_anno_2 <- column_anno_numerical[[a]]
      tmp_col_anno_2 <- tmp_col_anno_2[base::colnames(mat_heatmap)]
      if(base::is.null(anno_list)){
        anno_list <- ComplexHeatmap::HeatmapAnnotation(cont_anno = ComplexHeatmap::anno_boxplot(tmp_col_anno_2, height = grid::unit(1, "cm")),
                                                                   which = "column",
                                                                   annotation_name_side = "right",
                                                                   gap = grid::unit(2, "mm"),
                                                                   annotation_name_rot = 0,
                                                                   annotation_name_gp = grid::gpar(fontsize = 8),
                                                                   annotation_label = base::names(column_anno_numerical)[a], show_legend = F)
      }else{
        anno_list <-  ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::HeatmapAnnotation(cont_anno = ComplexHeatmap::anno_boxplot(tmp_col_anno_2, height = grid::unit(1, "cm")),
                                                                   which = "column",
                                                                   annotation_name_side = "right",
                                                                   gap = grid::unit(2, "mm"),
                                                                   annotation_name_rot = 0,
                                                                   annotation_name_gp = grid::gpar(fontsize = 8),
                                                                   annotation_label = base::names(column_anno_numerical)[a], show_legend = F), direction = "vertical")
      }
      

    }
  }


  all_conditions <- NULL

  # if(!is.null(column_anno_categorical) | !is.null(column_anno_numerical)){
    for(setnum in 1:base::length(hcobject[["layers"]])){
      all_conditions <- base::c(all_conditions, base::as.character(dplyr::pull(data[[base::paste0("set", setnum, "_anno")]], group)))
    }
    all_conditions <- base::table(all_conditions) %>%
      base::as.data.frame() %>%
      dplyr::filter(., all_conditions %in% base::colnames(mat_heatmap))
    all_conditions <- all_conditions[base::match(base::colnames(mat_heatmap), base::as.character(all_conditions$all_conditions)),]
    all_conditions <- base::paste0(all_conditions$all_conditions, "  [", all_conditions$Freq, "]")

    if(base::is.null(anno_list)){
        anno_list <- ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions))
    }else{
        anno_list <-  ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions)), direction = "vertical")
    }

    

  # }



  Cairo::Cairo(file = paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", file_name),
        width = 50,
        height = 30,
        pointsize=11,
        dpi=300,
        type = "pdf",
        units = "in")

  hm <- ComplexHeatmap::Heatmap(mat_heatmap,
                                right_annotation = ha,
                                col = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(base::length(base::seq(-2, 2, by = .1))),
                                clustering_distance_rows = "euclidean",
                                clustering_distance_columns = "euclidean",
                                clustering_method_rows = "complete",
                                clustering_method_columns = "complete",
                                cluster_columns = cluster_columns,
                                cluster_rows = cluster_rows,
                                column_names_rot = 90,
                                column_names_centered = F,
                                row_names_gp = grid::gpar(fontsize = 8),
                                column_names_gp = grid::gpar(fontsize = 8),
                                rect_gp = grid::gpar(col = "black"),
                                heatmap_legend_param = list(title = "", legend_height = grid::unit(3, "cm")), column_km = k)

  if(base::is.null(anno_list)){
    anno_list <- hm

  }else{
    anno_list <- ComplexHeatmap::add_heatmap(hm, anno_list)
  }

  hm_w_lgd <- ComplexHeatmap::draw(anno_list, annotation_legend_list = lgd_list, merge_legends = T,
                                   padding = grid::unit(c(2, 2, 2, 30), "mm"))

  grDevices::dev.off()

  print(hm_w_lgd)
  if(return_HM){
    return(hm_w_lgd)
  }

}