#' Plot Degree Distribution
#' 
#' A function that plot the degree distribution of each network after being cut with the set cutoff values.
#' @export

plot_deg_dist <- function(){
  for(x in 1:base::length(hcobject[["layers"]])){
    if(!hcobject[["cutoff_vec"]][x] %in% hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["cutoff_stats"]][["cutoff"]]){
      tmp_vec <- hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["cutoff_stats"]][["cutoff"]]
      hcobject[["cutoff_vec"]][x] <<- tmp_vec[base::which.min(base::abs(tmp_vec - hcobject[["cutoff_vec"]][x]))]
    }
    stats_calculated_optimal_cutoff <- hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["cutoff_stats"]][hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["cutoff_stats"]][["cutoff"]] == hcobject[["cutoff_vec"]][x], c("degree", "Probs")]
    
    stats <- hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["cutoff_calc_out"]][["cutoff_stats_concise"]]
    stats$cut_off <- base::rownames(stats)
    stats <- stats[stats[["cut_off"]] == hcobject[["cutoff_vec"]][x],]
    
    
    dd_plot_calculated_optimal <- ggplot2::ggplot(stats_calculated_optimal_cutoff, ggplot2::aes(x = base::log(degree), y = base::log(Probs))) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method="lm") +
      ggplot2::theme_bw() + 
      ggplot2::labs(x = "log(degree)",
                    y = "log(Probs)",
                    title = hcobject[["layers_names"]][x],
                    subtitle = base::paste0("Cut-off: ", hcobject[["cutoff_vec"]][x], "; R²: ", base::round(stats[1],3), "; no. edges: ",
                                            stats[2], "; no. nodes: ", stats[3], "; no. networks: ", stats[4]) )
      ggplot2::theme(plot.title = ggplot2::element_text(size = 14), plot.subtitle = ggplot2::element_text(size = 10)) 
    
    print(dd_plot_calculated_optimal)
    ggplot2::ggsave(base::paste0("Degree_distribution_plot_", hcobject[["layers_names"]][x], "_", hcobject[["cutoff_vec"]][x], ".pdf"),
          dd_plot_calculated_optimal, device = cairo_pdf, width = 10, height = 8, path = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]]))
    
  }
  
}
