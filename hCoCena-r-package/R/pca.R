#' Run And Plot PCA
#' 
#' Plots one PCA for each dataset.
#' @param which One of "all" (uses all genes), "topvar" (uses top most variant genes, cannot be used before running "Data processing part I" in the main markdown), 
#' 	or "network" (uses the genes present in network, cannot be run before finishing "Data processing part II" in the main markdown).
#' @param color_by Default is "none" (coloring all dots black), alternatively, set this to a vector of column names (c("name of col in annotation 1", "name of col in annotation 2", ...)), containing the groups by which you want to color.
#' @param ellipses A Boolean. Whether or not to add ellipses to the PCA plot. For details see documentation on factoextra::fviz_pca_ind. Default is FALSE.
#' @export


PCA <- function(which = "all", color_by = "none", ellipses = F){

  
  out <- base::lapply(1:base::length(hcobject[["layers"]]), function(x){
    
    if(which == "all"){
      res.pca <- stats::prcomp(base::t(hcobject[["data"]][[base::paste0("set", x, "_counts")]]), scale = F)
    }else if(which == "topvar"){
      res.pca <- stats::prcomp(base::t(hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["topvar"]]), scale = F)
    }else if(which == "network"){
      genes <- igraph::V(hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part2"]][["heatmap_out"]][["filt_cutoff_graph"]])$name
      dat <- hcobject[["data"]][[base::paste0("set", x, "_counts")]]
      dat <- dat[base::rownames(dat) %in% genes,]
      res.pca <- stats::prcomp(base::t(dat), scale = F)
    }else{
      message("Invalid parameter setting: which = ", which, " is not recognized.")
    }
    
    if(length(color_by) == 1){
      if(color_by == "none"){
        groups <- "none"
      }else{
        groups <- dplyr::pull(hcobject[["data"]][[base::paste0("set", x, "_anno")]], color_by)
      }
    }else{
      groups <- dplyr::pull(hcobject[["data"]][[base::paste0("set", x, "_anno")]], color_by[x])
    }
    
    if(base::length(base::unique(groups)) > base::length(ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(8))){
      my_palette <- grDevices::colorRampPalette(ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(8))(base::length(base::unique(groups)))
    }else(
      my_palette <- ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(base::length(base::unique(groups)))
    )
    
    g <- factoextra::fviz_pca_ind(res.pca,
                                   geom = "point",
                                   addEllipses = ellipses,
                                   habillage = groups,
                                   palette = my_palette,
                                   label = "none",
                                   title = base::paste0("PCA ", hcobject[["layers_names"]][x], " ", which),
                                   pointsize = 5, 
                                   invisible ="quali"
                                  )+
      ggplot2::scale_shape_manual(values = base::rep(19, base::length(base::unique(groups))))+
      ggplot2::theme_bw()
    
    graphics::plot(g)
    
    Cairo::CairoPDF(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/PCA_", which, "_", hcobject[["layers_names"]][x], ".pdf"), width = 10, height = 7)
    graphics::plot(g)
    grDevices::dev.off()

    return(g)
    
  })
  

  hcobject[["satellite_outputs"]][["pca"]] <<- out
}
