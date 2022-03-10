#' Plot Sample Distributions
#' 
#' The distribution of count values per sample for all datasets are visualized either as boxplots or frequency distributions.
#' @param plot_type Either "boxplot" or "freqdist", to create boxplots or frequency distributions, respectively.
#' @param log_2 A Boolean. Whether or not the counts should be logged before plotting.
#' @param plot A Boolean. Whether or not to to print the plots into the markdown in addition to saving them to the save folder. 
#'  If the total number of samples is very large, setting this parameter to FALSE is advised in order to prevent the knitted R Markdown from becoming too long.
#' @export

plot_sample_distributions <- function(plot_type = "boxplot", 
										log_2 = T, 
										plot = T){ #plot_type: one of "boxplot" or "freqdist"
  
  if(!"sample_distribution_plots" %in% base::list.dirs(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]],"/"))){
    
    base::dir.create(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]],"/","sample_distribution_plots"))
    
    }
  
  data <- list()
  
  for (i in 1:base::length(hcobject[["layers"]])){
    
    data[[i]] <- hcobject[["data"]][[base::paste0("set", i, "_counts")]]
    
  }
  
  # boxplot or frequency distibution?
  if(plot_type == "boxplot"){
    
    plt <- boxplot_from_list(data = data, log_2 = log_2, bool_plot = plot)
    
    if(plot == T){
      
      plot_list_of_plots(plt)
      
    }
    
  }
  else if(plot_type == "freqdist"){
    
    plt <- freqdist_plot_from_list(data = data, log_2 = log_2, bool_plot = plot)
    
    if(plot == T){
      
      plot_list_of_plots(plt)
      
    }
    
  }else{
    
    print("plot_type should be either 'boxplot' or 'freqdist'.")
    
  }

  hcobject[["satellite_outputs"]][["sample_distribution_plots"]] <<- plt  

}

