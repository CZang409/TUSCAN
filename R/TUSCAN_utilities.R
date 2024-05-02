#####################################################################
# Package: TUSCAN
# Date : 2024-04-24
# Title : Tumor Segmentation and Classification Analysis in Spatial Transcriptomics
# Authors: Chenxuan Zang
# Contacts: czang@mdanderson.org
#          University of Texas MD Anderson, Department of Biostatistics
######################################################################

#' The TUSCAN Class
#'
#'Slots in the TUSCAN object include:
#'
#' @slot counts.data The raw expression count matrix. Rows are genes, columns are spots.
#' @slot cnv.data The copy number matrix inferred from the raw count matrix. Rows are genes, columns are spots.
#' @slot project Name of the project.
#' @slot location The spatial coordinates of the spots, with two columns representing the \eqn{x} and \eqn{y} coordinates of the spatial location.
#' @slot image The histology image information, including the pixel location of each spot on the original histology image, and color information.
#' @slot autoencoder The outputs from autoencoder.
#' @slot pre_cluster  The preliminary clustering results that are used to select a reference cluster later.
#' @slot tumor_cluster Tumor prediction, indicating each spot as 'tumor' or 'normal'.
#' @slot annotation Annotations of all spots, with two columns representing the spot barcode and the spot clustering name. If the user do not provide an annotation file, the \code{pre_cluster} result will be used as the annotations.
#' @slot gene_order The gene order file, contains chromosome, start, and stop position for each gene.
#' @slot reference_grouped_spot_indices The indices of the reference spots.
#' @slot observation_grouped_spot_indices The indices of the other spots.
#'
#' @return
#' @import magick
#' @export
#'
#' @examples
#'

setClass("TUSCAN", slots=list(
  counts.data = "ANY",
  cnv.data = "ANY",
  project = "character",
  location = "data.frame",
  image = "ANY",
  autoencoder = "ANY",
  pre_cluster = "ANY",
  tumor_cluster = "ANY",
  annotation = "ANY",
  gene_order= "data.frame",
  reference_grouped_spot_indices = "list",
  observation_grouped_spot_indices = "list"

))


#' Create a TUSCAN object
#'
#' Create a TUSCAN object.
#' @param counts Raw spatial transcriptomics gene expression count matrix, the dimension is \eqn{m \times n}, where \eqn{m} is the number of genes and \eqn{n} is the number of spots.
#' @param location Dataframe containing the spatial coordinates of all spots, where each row represents a spot and two columns represent the \eqn{x} and \eqn{y} coordinates of the spatial location, respectively.
#' @param project Name of the project.
#' @param img Histology image input.
#' @param r Image smoothing parameter: draw a square with a side length of (\eqn{2r+1}) centered on each spot. Apply the average color within the square as the new RGB values. The default value is 49.
#'
#' @return Return a preliminary TUSCAN object.
#' @import magick
#' @export
#'
#' @examples

createTUSCANObject <- function(counts,location,project="TUSCAN",img, r=49){

  ## check dimension
  if(ncol(counts)!=nrow(location)){
    stop("The number of spots in counts and location should be consistent.")
  }

  ## check data order
  if(!identical(colnames(counts), rownames(location))){
    stop("The column names of counts and row names of location should be should be matched.")
  }

  pixel_data <- image_data(img)
  # convert rgb values to 0-255
  image_values <- as.integer(pixel_data)

  # create an empty dataframe
  image_df <- data.frame(barcode = numeric(nrow(location)),
                         x = numeric(nrow(location)),
                         y = numeric(nrow(location)),
                         R = numeric(nrow(location)),
                         G = numeric(nrow(location)),
                         B = numeric(nrow(location)),
                         RS = numeric(nrow(location)),
                         GS = numeric(nrow(location)),
                         BS = numeric(nrow(location))
  )

  max_x = max(location$x)
  max_y =max(location$y)
  radius = r

  for(i in 1:nrow(location)){
    image_df$barcode[i] = row.names(location)[i]
    image_df$y[i] = location$y[i]
    image_df$x[i] = location$x[i]
    # original rgb values
    image_df$R[i] = image_values[location$y[i],location$x[i],1]
    image_df$G[i] = image_values[location$y[i],location$x[i],2]
    image_df$B[i] = image_values[location$y[i],location$x[i],3]
    # take average rgb values around the spot
    image_df$RS[i] <- round(mean(image_values[max(0,location$y[i]-radius):min(max_y,location$y[i]+radius),
                                              max(0,location$x[i]-radius):min(max_x,location$x[i]+radius),
                                              1]))
    image_df$GS[i] <- round(mean(image_values[max(0,location$y[i]-radius):min(max_y,location$y[i]+radius),
                                              max(0,location$x[i]-radius):min(max_x,location$x[i]+radius),
                                              2]))
    image_df$BS[i] <- round(mean(image_values[max(0,location$y[i]-radius):min(max_y,location$y[i]+radius),
                                              max(0,location$x[i]-radius):min(max_x,location$x[i]+radius),
                                              3]))
  }

  image_df$gray <- image_df$RS*0.299+
    image_df$GS*0.587+
    image_df$BS*0.114

  object <- new(

    Class = "TUSCAN",
    counts.data = counts,
    project = project,
    location = location,
    image = image_df

  )

  return(object)
}
