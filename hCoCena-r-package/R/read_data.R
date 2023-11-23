
#' Function To Read All Count And Annotation Files
#' 
#' The function loads the count and annotation data for each layer and saves it in the hCoCena-Object's "data" slot.
#' @param sep_counts The separator of the count files. Default is tab separated files. Ignore when loading data from objects instead of files.
#' @param sep_anno The separator of the annotation files. Default is tab separated files. Ignore when loading data from objects instead of files.
#' @param gene_symbol_col A String. Name of the column that contains the gene symbols. Ignore when loading data from objects instead of files.
#' @param sample_col A String. Name of the column that contains the sample IDs. Ignore when loading data from objects instead of files.
#' @param count_has_rn A Boolean. Whether or not the count file has rownames. Default is TRUE. Ignore when loading data from objects instead of files.
#' @param anno_has_rn A Boolean. Whether or not the annotation file has rownames. Default is TRUE. Ignore when loading data from objects instead of files.
#' @export

read_data <- function(sep_counts = "\t", 
						sep_anno = "\t", 
						gene_symbol_col = NULL, 
						sample_col = NULL, 
						count_has_rn = TRUE, 
						anno_has_rn = TRUE){

	base::options(dplyr.summarise.inform = F)
  
	data <- list()

	# iterate over layers:
	for (x in 1:base::length(hcobject[["layers"]])){
		# read expression data:

		# check if provided string refers to an existing object, if so load the object instead of reading from file:
		if(base::exists(hcobject[["layers"]][[base::paste0("set",x)]][1])){
		  
		  data[[base::paste0("set", x, "_counts")]] <- base::get(hcobject[["layers"]][[base::paste0("set",x)]][1])
		  
		}else{
		  
		  if(base::is.null(gene_symbol_col)){
		  	stop("You must provide the 'gene_symbol_col' parameter.")
		  }
		  data[[base::paste0("set", x, "_counts")]] <- read_expression_data(file = base::paste0(hcobject[["working_directory"]][["dir_count_data"]], hcobject[["layers"]][[base::paste0("set",x)]][1]), rown = count_has_rn,
		                                                              sep = sep_counts, gene_symbol_col = gene_symbol_col)
		}
	  # check for zero-variance genes
	    var.df <- rank_variance(data[[base::paste0("set", x, "_counts")]])
	    if(any(var.df$variance == 0)){
	    message(base::paste0("Detected genes with 0 variance in dataset ", x, "."))
	    data[[base::paste0("set", x, "_counts_unfiltered")]] <- data[[base::paste0("set", x, "_counts")]]
	    data[[base::paste0("set", x, "_counts")]] <- data[[base::paste0("set", x, "_counts")]] %>%
	      dplyr::filter(!(row.names(.) %in% (var.df %>% dplyr::filter(variance == 0) %>% dplyr::pull(gene))))
	    message(base::paste0((nrow(data[[base::paste0("set", x, "_counts_unfiltered")]])-nrow(data[[base::paste0("set", x, "_counts")]])), " gene(s) were removed from dataset ", x, "."))
	    } 
	  

		# read annotation data:

		# check if provided string refers to an existing object, if so load the object instead of reading from file:
		if(base::exists(hcobject[["layers"]][[base::paste0("set",x)]][2])){
		  anno <- base::get(hcobject[["layers"]][[base::paste0("set",x)]][2])
		  
		  anno[] <- base::lapply(anno, base::factor)
		  
		  data[[base::paste0("set", x, "_anno")]] <- anno
		  
		}else{
		  if(base::is.null(sample_col)){
		  	stop("You must provide the 'sample_col' parameter.")
		  }
		  data[[base::paste0("set", x, "_anno")]] <- read_anno(file = base::paste0(hcobject[["working_directory"]][["dir_annotation"]], hcobject[["layers"]][[base::paste0("set",x)]][2]), rown = anno_has_rn,
		                                                 sep = sep_anno, sample_col = sample_col)
		  
		}

		# check if samples match between annotation and counts:
		if(!base::ncol(data[[base::paste0("set", x, "_counts")]]) == base::nrow(data[[base::paste0("set", x, "_anno")]])){
		  stop(base::paste0("The count table has ",  base::ncol(data[[base::paste0("set", x, "_counts")]]), " columns but the annotation has ", 
		              base::nrow(data[[base::paste0("set", x, "_anno")]]), " rows. These values are required to be the same since they
		                 should correspond to the number of samples. THE LOADING OF THE DATA WILL BE TERMINATED."))
		}else if(!base::all(base::as.character(base::colnames(data[[base::paste0("set", x, "_counts")]]))  %in% base::as.character(base::rownames(data[[base::paste0("set", x, "_anno")]])))){
		  stop(base::paste0("The column names of the count file do not all match the rownames of the annotation. Please make sure they contain the same samples."))
		}else{
		  data[[base::paste0("set", x, "_anno")]] <- data[[base::paste0("set", x, "_anno")]][base::colnames(data[[base::paste0("set", x, "_counts")]]),]
		}

	}

	hcobject[["data"]] <<- data
}
