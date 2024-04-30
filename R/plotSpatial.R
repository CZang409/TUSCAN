
#' Spatial plotting functions
#'
#' @param obj TUSCAN object.
#' @param cluster Options: "cluster","tumor".
#' @param color Optional color palette can be customized by user.
#' @param size Point size.
#'
#' @return Return a spatial plot
#' @import ggplot2
#' @export
#'
#' @examples
#'

plotSpatial <- function(obj, class = c("cluster", "tumor"), color = NULL, size){

  class <- match.arg(class)

  if (is.null(color)) {
    color <- c(RColorBrewer::brewer.pal(8, "Set1"),RColorBrewer::brewer.pal(8, "Set2"))
  }

  df <- obj@image
  df[[class]] <- as.factor(df[[class]])

  p <- ggplot(df, aes_string(x = "x", y = "y", color = class)) +
    geom_point(size = size) +
    scale_y_reverse()+
    scale_color_manual(values = color) +
    labs(color="Cluster") +
    coord_equal() +
    theme_void()

  return(p)

}
