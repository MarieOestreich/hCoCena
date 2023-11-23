#' Meta Data Plot
#' 
#' Visualize the meta data with respect to a grouping of choice, e.g., how outcome or age are distributed across the groups of your variable of interest. 
#'  In case of categorical data, stacked bar plots will be plotted, in case of numerical data it will be box plots.
#' @param set An integer. The number of the dataset for which the plot should be  created.
#' @param group_col A string. The name of the column in the annotation of datset 'set' that contains your groups (e.g. conditions) that are to be plotted along the x-axis.
#' @param meta_col A string. The name of the column in the annotation of datset 'set' that contains the information that is to be inspected (e.g. outcome or age).
#' @param type Either "cat", if "meta_col" is categorical, or "num", if "meta_col" is numerical.
#' @param cols User-defined color vector.
#' @export


meta_plot <- function(set, group_col = NULL, meta_col = NULL, type = "cat", cols = NULL){

  anno <- hcobject[["data"]][[base::paste0("set", set, "_anno")]]

  if(base::is.null(group_col)){
    stop("No annotation data columns for group reference provided (parameter 'group_col').")
  }
  if(base::is.null(meta_col)){
    stop("No annotation data columns for plotting provided (parameter 'meta_col').")
  }
  
  
  if(!group_col %in% base::colnames(anno)){
    stop(base::paste0("*", group_col, "* is not a valid column name of the annotation table. Check for spelling errors and note that this parameter is case sensitive."))
  }
  if(!meta_col %in% base::colnames(anno)){
    stop(base::paste0("*", meta_col, "* is not a valid column name of the annotation table. Check for spelling errors and note that this parameter is case sensitive."))
  }
  
  plot_df <- NULL
  
  if(type == "cat"){
    base::lapply(base::unique(dplyr::pull(anno, group_col)), function(x){
      tmp <- anno[anno[,group_col] == x, ]
      total <- base::nrow(tmp)
      tb <- base::table(dplyr::pull(tmp, meta_col)) %>% base::as.data.frame()
      tb$Freq <- tb$Freq/base::sum(tb$Freq)
      plot_df <<- base::rbind(plot_df, base::data.frame(group = base::paste0(x, " [", total, "]"), meta = tb$Var1, value = tb$Freq))
    })
    

    if(base::is.null(cols)){
      if(base::length(base::unique(plot_df$meta)) > base::length(ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(8))){
        my_palette <- grDevices::colorRampPalette(ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(8))(base::length(base::unique(plot_df$meta)))
      }else(
        my_palette <- ggsci::pal_nejm(palette = c("default"), alpha = 1)(base::length(base::unique(plot_df$meta)))
      )
    }else{
      my_palette <- cols
    }
    
    p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = group, y = value, fill = meta)) +
      ggplot2::geom_bar(stat = "identity")+ 
      ggplot2::scale_fill_manual(values = my_palette)+
      ggplot2::theme_bw()
    graphics::plot(p)
  }
  
  if(type == "num"){
    base::lapply(base::unique(dplyr::pull(anno, group_col)), function(x){
      tmp <- anno[anno[, group_col] == x, ]
      total <- base::nrow(tmp)
      plot_df <<- base::rbind(plot_df, base::data.frame(group = base::paste0(x, " [", total, "]"), value = base::as.numeric(base::as.character(dplyr::pull(tmp, meta_col)))))
    })
    
    if(base::is.null(cols)){
      if(base::length(base::unique(plot_df$group)) > base::length(ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(8))){
        my_palette <- grDevices::colorRampPalette(ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(8))(base::length(base::unique(plot_df$group)))
      }else{
        my_palette <- ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(base::length(base::unique(plot_df$group)))
      }
    } else {
      my_palette <- cols
    }
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = value)) + 
      ggplot2::geom_boxplot(ggplot2::aes(fill = group))+ 
      ggplot2::scale_fill_manual(values = my_palette)+
      ggplot2::theme_bw()

    graphics::plot(p)
  }

}