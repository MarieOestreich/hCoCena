
#' Plot Cluster Heatmap
#' 
#' Plots a heatmap with sample groups as columns and gene clusters as rows. The cells are coloured according to the mean GFC of a given cluster in the respective sample group.
#'  If categorical or numerical metadata annotations or user-defined enrichments have been created with satellite functions, they will be incorporated into the heatmap as column/row annotations.
#'  For the available options check out the satellite functions.
#' @param col_order Defines the order in which the sample groups (conditions) appear in the heatmap. 
#'  Accepts a vector of strings giving the conditions in their desired order, default is NULL (automatic order, determined by clustering).
#'  If the parameter cluster_columns is not set to FALSE, this order will be overwritten.
#' @param row_order Like col_order but with cluster names.
#' @param cluster_columns A Boolean, whether or not to cluster the columns of the heatmap. Default is TRUE, overwrites col_order.
#' @param cluster_rows Like cluster_columns but for rows.
#' @param k The resulting cluster_columns tree is cut into k groups. Default is 0 (no cutting).
#' @param return_HM A Boolean whether of not to return the ComplexHeatmap object to hcobject$integrated_output$cluster_calc$heatmap_cluster in addition to plotting it. Default is TRUE.
#' @param cat_as_bp  A vector of Booleans which length is equivalent to the number of categorical meta data annotations created with the satellite functions. 
#'  Each Boolean states whether or not the categorical variable should be annotated as a bar plot (TRUE) or as a line plot (FALSE). 
#'  If you did not perform any meta data annotation, ignore this parameter. 
#' @param file_name A string giving the name of the file (with .pdf ending) to which the heatmap should be written. Default is "module_heatmap.pdf".
#' @export



plot_cluster_heatmap <- function(col_order = NULL, 
                                       row_order = NULL, 
                                       cluster_columns = T,
                                       cluster_rows = T, 
                                       k = 0, 
                                       return_HM = T, 
                                       cat_as_bp = NULL, 
                                       file_name = "Heatmap_modules.pdf"){

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
  
  # get categorical and numerical sample group annotations if they exist:
  if("column_annos_categorical" %in% base::names(hcobject[["satellite_outputs"]])){
    column_anno_categorical <- hcobject[["satellite_outputs"]][["column_annos_categorical"]]
  }else{
    column_anno_categorical <- NULL
  }
  
  if("column_annos_numerical" %in% base::names(hcobject[["satellite_outputs"]])){
    column_anno_numerical <- hcobject[["satellite_outputs"]][["column_annos_numerical"]]
  }else{
    column_anno_numerical <- NULL
  }
  
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
      c_GFCs <- dplyr::filter(hcobject[["integrated_output"]][["GFC_all_layers"]], Gene %in% genes)
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
      c_GFCs <- dplyr::filter(hcobject[["integrated_output"]][["GFC_all_layers"]], Gene %in% genes)

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


  base::colnames(mat_heatmap) <- base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])[1:(base::ncol(hcobject[["integrated_output"]][["GFC_all_layers"]])-1)]

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
      tmp_colour <- grDevices::colorRampPalette(c("#332288", "#117733", "#44aa99", "#88ccee", "#cc6677", "#aa4499", "#882255"))(base::ncol(column_anno_categorical[[a]]))
      
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
                                                                     annotation_label = base::names(column_anno_categorical)[a], direction = c("vertical"))
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
                                                                     annotation_label = base::names(column_anno_categorical)[a]), direction = c("vertical"))
        }
        
        
        
      }else{
        if(base::is.null(anno_list)){
          anno_list <-ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_lines(column_anno_categorical[[a]] %>% base::as.matrix(),width = grid::unit(2, "cm"),
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
                                                                     annotation_label = base::names(column_anno_categorical)[a]), direction = c("vertical"))
        }
        
      }
      
      

      lgd_list <- rlist::list.append(lgd_list, ComplexHeatmap::Legend(labels = base::colnames(column_anno_categorical[[a]]%>%base::as.matrix()), 
                                                                      title = base::names(column_anno_categorical)[a],
                                                                      legend_gp = grid::gpar(col = tmp_colour),
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
                                                                   annotation_label = base::names(column_anno_numerical)[a], show_legend = F), direction = c("vertical"))
      }
      

    }
  }


  all_conditions <- NULL

  # if(!is.null(column_anno_categorical) | !is.null(column_anno_numerical)){
    for(setnum in 1:base::length(hcobject[["layers"]])){
      all_conditions <- base::c(all_conditions, base::as.character(dplyr::pull(hcobject[["data"]][[base::paste0("set", setnum, "_anno")]], hcobject[["global_settings"]][["voi"]])))
    }
    all_conditions <- base::table(all_conditions) %>%
      base::as.data.frame() %>%
      dplyr::filter(., all_conditions %in% base::colnames(mat_heatmap))
    all_conditions <- all_conditions[base::match(base::colnames(mat_heatmap), base::as.character(all_conditions$all_conditions)),]
    all_conditions <- base::paste0(all_conditions$all_conditions, "  [", all_conditions$Freq, "]")
    if(base::is.null(anno_list)){
      anno_list <- ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions))

    }else{
      anno_list <-  ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions)), direction = c("vertical"))
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
                                #col = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(base::length(base::seq(-2, 2, by = .1))),
                                col = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(51),
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
                                heatmap_legend_param = list(title = "GFC", legend_height = grid::unit(3, "cm")), column_km = k)
  if(base::is.null(anno_list)){
    anno_list <- hm 
  }else{
    anno_list <- ComplexHeatmap::add_heatmap(hm, anno_list, direction = c("vertical"))
  }

  hm_w_lgd <- ComplexHeatmap::draw(anno_list, annotation_legend_list = lgd_list, merge_legends = T,
                                   padding = grid::unit(c(2, 2, 2, 30), "mm"))

  grDevices::dev.off()

  print(hm_w_lgd)
  if(return_HM){
    hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]] <<- hm_w_lgd
  }

}
