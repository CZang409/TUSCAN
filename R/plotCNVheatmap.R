#' Create a color gradient palette
#'
#' Create a color gradient palette.
#' @return
#'
#' @keywords internal
#' @noRd
#'

get_group_color_palette <- function() {
  return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")))
}


#' Plot the copy number matrix as a heatmap, with spots as rows and genes as columns, ordered according to chromosome
#'
#' Visualize the CNV profile by heatmap.
#' @param object TUSCAN object.
#' @param color_belt Optional color palette for plotting.
#' @param title Plot title.
#'
#' @return Create a PNG or PDF image file.
#' @import ComplexHeatmap
#' @importFrom grid gpar
#' @export
#'
#' @examples
plotCNVheatmap <- function(object,color_belt=NULL,output_format=c("png", "pdf"),title="heatmap"){

  output_format <- match.arg(output_format)

  col_fun = circlize::colorRamp2(c(mean(object@cnv.data)-3*sd(object@cnv.data), 1, mean(object@cnv.data)+3*sd(object@cnv.data)), c("blue", "white", "red"))
  expr <- object@cnv.data

  # normal
  normal_loc <- object@reference_grouped_spot_indices
  ref_cluster = names(normal_loc)
  #color row
  row_color_ref = data.frame(cluster=rep(ref_cluster[1],length(normal_loc[[1]])))
  if(length(ref_cluster)>1){
    for(k in 2:length(ref_cluster)){
      tmp = data.frame(cluster=rep(ref_cluster[k],length(normal_loc[[k]])))
      row_color_ref = rbind(row_color_ref, tmp)
    }
  }
  normal_loc <- unlist(normal_loc)

  # other
  other_loc <- object@observation_grouped_spot_indices
  obs_cluster = names(other_loc)
  ##row color
  row_color_obs = data.frame(cluster=rep(obs_cluster[1],length(other_loc[[1]])))
  if(length(obs_cluster)>1){
    for(k in 2:length(obs_cluster)){
      tmp = data.frame(cluster=rep(obs_cluster[k],length(other_loc[[k]])))
      row_color_obs = rbind(row_color_obs, tmp)
    }
  }
  other_loc <- unlist(other_loc)

  sub_geneFile <-  object@gene_order
  #identical with expr
  identical(rownames(expr),rownames(sub_geneFile)) # TRUE

  # expr
  norm_expr <- as.data.frame(expr[,normal_loc])
  norm_expr$chr <- as.factor(sub_geneFile[,'chr'])
  other_expr <- as.data.frame(expr[,other_loc])
  other_expr$chr <- as.factor(sub_geneFile[,'chr'])

  new_cluster = other_expr$chr
  if(is.null(color_belt)){
    color <- get_group_color_palette()(length(unique(other_expr$chr)))
  }else{
    color = color_belt[1:length(unique(other_expr$chr))]
  }
  # annotation
  top_color <- HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = color),
                         labels = levels(new_cluster),
                         labels_gp = gpar(cex = 1,
                                          col = "black"
                         )))
  #row color
  all_cluster = unique(c(ref_cluster,obs_cluster))

  if(is.null(color_belt)){
    rowcolor = get_group_color_palette()(length(all_cluster))
  }else{
    rowcolor = color_belt[1:length(all_cluster)]
  }
  names(rowcolor) = all_cluster

  row_color1 = rowAnnotation(
    cluster = row_color_ref[,1],
    col = list(cluster = rowcolor[ref_cluster]),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  row_color2 = rowAnnotation(
    cluster = row_color_obs[,1],
    col = list(cluster = rowcolor[obs_cluster]),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  lgd = Legend(col_fun = col_fun, title = "Expression Distribution",
               legend_height = unit(6, "cm"),
               grid_width = unit(2, "cm"),
               title_position = "leftcenter-rot",
               title_gp = gpar(fontsize = 20),
               labels_gp = gpar(fontsize = 15))

  lgd_normal = Legend(labels = ref_cluster,
                      legend_gp = gpar(fill = rowcolor[ref_cluster]),
                      title = "Normal",
                      nrow = 1,
                      legend_height = unit(1, "cm"),
                      grid_width = unit(1, "cm"),
                      title_gp = gpar(fontsize = 18),
                      labels_gp = gpar(fontsize = 16),
                      gap = unit(6, "mm")
  )

  lgd_other = Legend(labels = obs_cluster,
                     legend_gp = gpar(fill = rowcolor[obs_cluster]),
                     title = "Others",
                     nrow = 1,
                     legend_height = unit(1, "cm"),
                     grid_width = unit(1, "cm"),
                     title_gp = gpar(fontsize = 18),
                     labels_gp = gpar(fontsize = 16),
                     gap = unit(6, "mm")
  )

  pd = packLegend(lgd_normal, lgd_other,
                  row_gap = unit(6, "mm"))

  n <- t(other_expr[,-ncol(other_expr)])
  m <- t(norm_expr[,-ncol(norm_expr)])

  ht_normal = Heatmap(as.matrix(m),
                      col = col_fun,
                      heatmap_legend_param = list(
                        col_fun = col_fun, title = "Expression Distribution",
                        legend_height = unit(6, "cm"),
                        grid_width = unit(2, "cm"),
                        title_position = "leftcenter-rot",
                        title_gp = gpar(fontsize = 20),
                        labels_gp = gpar(fontsize = 15),
                        padding = unit(10, "cm")
                      ),
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_column_names = FALSE,
                      show_row_names = FALSE,
                      column_split = new_cluster,
                      width = unit(50, "cm"),
                      height = unit(10, "cm"),
                      row_title = "Normal Control",
                      row_title_side = c("left"),
                      row_title_rot = 90,
                      left_annotation = row_color1,
                      row_title_gp = gpar(fontsize = 25),
                      column_title = "CNV Heatmap",
                      column_title_gp = gpar(fontsize = 35, fontface = "bold"),
                      row_split = row_color_ref[,1],
                      border = TRUE,
                      row_gap = unit(0, "mm"),
                      column_gap = unit(0, "mm")
  )

  ht_other = Heatmap(as.matrix(n),
                     col = col_fun,
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     show_column_names = FALSE,
                     show_row_names = FALSE,
                     column_split = new_cluster,
                     width = unit(50, "cm"),
                     height = unit(30, "cm"),
                     show_heatmap_legend=FALSE,
                     top_annotation = top_color,
                     left_annotation = row_color2,
                     row_title = "Other Spots",
                     row_title_side = c("left"),
                     row_title_rot = 90,
                     row_title_gp = gpar(fontsize = 25),
                     column_title = "I am a column title at the bottom",
                     column_title_side = "bottom",
                     row_split = row_color_obs[,1],
                     border = TRUE,
                     row_gap = unit(0, "mm"),
                     column_gap = unit(0, "mm")
  )

  # plot

  if(output_format=="png"){

    png(paste(title,".png",sep=""),width = 25/7*480,height = 20/7*480)

    draw(ht_normal %v% ht_other,
        annotation_legend_list = pd,
        annotation_legend_side = "bottom",
        heatmap_legend_side = "right"
    )

    dev.off()

  }else if(output_format=="pdf"){

    pdf(paste(title,".pdf",sep=""),width = 25,height = 20)

    draw(ht_normal %v% ht_other,
         annotation_legend_list = pd,
         annotation_legend_side = "bottom",
         heatmap_legend_side = "right"
    )

    dev.off()
  }

}
