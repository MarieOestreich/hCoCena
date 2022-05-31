#' Run First Part Of The Gene Expression Analysis
#' 
#' This function executes the frist part of the data processing procedure. It leads up to choosing the correlation cut-off for each layer.
#' 	All datasets will be filtered for their most variant genes as defined in the layer-specific settings. 
#' 	After this filtering step, the pair-wise Pearson correlation coefficients for all pairs of genes are calculated. 
#' 	Correlations that are negative or that have an associated p-value higher than 0.05 are immediately discarded.
#' 	Next, a set of statistics will be calculated for the set range of cut-off values that aim to facilitate the cut-off choice. 
#' 	This includes determining the number of graph components resulting from creating a network when cutting the data with the respective cut-off, 
#' 	as well as the number of nodes and edges this network comprises. 
#' 	The last parameter that is evaluated is the R²-value of the data to a linear regression through the logged degree distribution for the given network.
#' @param padj A String. Defines the method to be used for p-value adjustment. Valid values are "none" (default) or values for "method" in stats::p.adjust.
#' @param export A Boolean. If TRUE, correlation values and p-values will be exported. 
#' 	This can save time if you plan on re-running the analysis since computing pari-wise correlations is a bottleneck of the analysis. Default is FALSE.
#' @param import A list. Each slot in the list corresponds to one of the layers (datasets) and is either a vector of two strings (1. path to file holding the correlation matrix and 
#'  2. Path to the file holding the p-value matrix) or NA. A list slot is set to NA if for that layer you do not want to import a pre-calculated correlation matrix. 
#'  The files do not necessarily have to be exported from a previous run, but can have any kind of origin (created with a different program or method). 
#'  For compatibility it is only important, that it is a whitespace-separated text (.txt) file containing a symmetric, numeric matrix. 
#'  The first line has to be gene names, there must be no row names, since the first line will be used for column names and row names. 
#'  Also, you must provide a matrix with correlation values AND a matrix with corresponding p-values, where cells in the matrices correspond to each other (only a correlation matrix will not be sufficient).
#'  Default is NULL.
#' @param bayes Sánchez-Taltavull et al. (2016) suggest superiority of Bayesian correlation analysis to Pearson correlation in some cases. 
#' 	Therefore, the Pearson correlation values can be weighted with Bayesian correlation values. To do so, set the “bayes”-parameter to TRUE. Default is FALSE, using only Pearson correlations.
#' @param alpha A numeric value from 0 to 1. Allows to adjust the strength of the Bayes weighting: For alpha = 0 the Pearson correlation values remain unaltered, for alpha = 1 the Pearson correlation value and the Bayesian correlation value contribute equally to the final correlation.
#' @param prior An integer, either 2 or 3, using prior 2 or 3 for the Bayes weighting as described in "Bayesian correlation analysis for sequence count data" by Sánchez-Taltavull et al. (2016).
#' @param corr_method Default is "pearson", but can alternatively be set to "spearman" or "rho" in case of single-cell data (according to Skinnider et al., https://www.nature.com/articles/s41592-019-0372-4).
#' @export

run_expression_analysis_1 <- function(padj = "none",
										export = F, 
										import = NULL, 
										bayes = F, 
										prior = 2, 
										alpha = 0.5,
										corr_method = "pearson"){
	if(!corr_method %in% c("pearson","spearman", "rho")){
		stop("Parameter 'corr_method' must be either 'pearson', 'spearman' or 'rho'.")
	}
	# run first part of expression analysis for each data layer:
	for(x in 1:base::length(hcobject[["layers"]])){
		hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]] <<- run_expression_analysis_1_body(x, 
		                                                                                        bayes = bayes, 
		                                                                                        prior = prior, 
		                                                                                        alpha = alpha, 
		                                                                                        padj = padj,
		                                                                                        export = export,
		                                                                                        import = import,
		                                                                                        corr_method = corr_method)
	}

}

