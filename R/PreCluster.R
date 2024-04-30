
#' Perform dimension reduction by autoencoder
#'
#' Perform dimension reduction by autoencoder.
#' @param object TUSCAN object.
#' @param nhvg Number of highly variable genes to run autoencoder.
#' @param d Number of dimensions compressed by autoencoder.
#' @param epochs Number of training iterations.
#' @param batch_size Number of training samples utilized during one iteration.
#' @param seed Optional numerical value. Set a random seed for reproducible results.
#'
#' @return Return a processed TUSCAN object adding low-dimensional components.
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures VariableFeatures
#' @export
#'
#' @examples

reduceDim <- function(object,nhvg=2000,d=20,epochs=50,batch_size=64,seed=NULL){

  Seurat_obj <- CreateSeuratObject(counts = object@counts.data, min.cells = 0,
                                   min.genes = 0, project = "TUSCAN")
  Seurat_obj <- NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  hvgs <- FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
  # Identify the 2000 most highly variable genes (hvgs)
  top2k <- VariableFeatures(hvgs)

  ####### Pre-clustering #########
  ### clustering based on autoencoder + scaled RGB

  log_mat <- Matrix::t(Seurat_obj@assays$RNA@data)
  hvg_mat <- log_mat[,top2k]

  # autoencoder
  print("Start dimension reduction by autoencoder...")

  features_matrix <- Matrix::as.matrix(hvg_mat)
  input_dim <- nhvg
  encoding_dim <- d # dim after dimension decrease

  set.seed(seed)
  tensorflow::set_random_seed(seed)
  #tensorflow::tf$random$set_seed(seed)
  input_layer <- keras::layer_input(input_dim)

  # set different number of hidden layers
  encoder_layer <- input_layer %>%
    keras::layer_dense(units = 1500, activation = "relu") %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dense(units = 1000, activation = "relu") %>%
    keras::layer_dense(units = 500, activation = "relu") %>%
    keras::layer_dense(units = 100, activation = "relu") %>%
    keras::layer_dense(units = encoding_dim) # 2 dimensions for the output layer

  decoder_layer <- encoder_layer %>%
    keras::layer_dense(units = 100, activation = "relu") %>%
    keras::layer_dense(units = 500, activation = "relu") %>%
    keras::layer_dense(units = 1000, activation = "relu") %>%
    keras::layer_dense(units = 1500, activation = "relu") %>%
    keras::layer_dense(units = input_dim)

  autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder_layer)

  autoencoder_model %>% keras::compile(
      loss='mean_squared_error',
      optimizer='adam',
      metrics = c('accuracy')
    )

  #set.seed(seed)
  autoencoder_model %>%
    keras::fit(features_matrix,
               features_matrix,
               epochs=epochs,
               batch_size = batch_size,
               shuffle=TRUE,
               verbose=0
    )

  encoder <- keras::keras_model(inputs = autoencoder_model$input, outputs = encoder_layer)
  encoded_data <- predict(encoder, features_matrix)
  object@autoencoder <- encoded_data
  print("Dimension reduction finished...")

  return(object)

}

#' Rescale a numerical vector to a new range
#'
#' This function performs a linear transformation on a numerical vector \code{x}, adjusting its range to \code{(gmin,gmax)}.
#' @param x A numeric vector that will be rescaled to a new specified range.
#' @param gmin Minimum value of the new range.
#' @param gmax Maximum value of the new range.
#'
#' @return Return a rescaled vector with its range adjusted to \code{(gmin,gmax)}.
#'
#' @keywords internal
#' @noRd
#'

rescale_to_range <- function(x,gmin,gmax){
  rescaled_x <- (x-min(x))/(max(x)-min(x))*(gmax-gmin)+gmin
  return(rescaled_x)
}

#' Preliminary clustering integrating spatial gene expression data and histology image
#'
#' This function performs a preliminary clustering of all spots based on spatial gene expression matrix and histology image.
#' The aim of this step is to select a subset of normal spots as the reference to infer copy number profile later.
#' We recommend using a larger number of clusters to make the selected reference group is purer.
#' We recommend that the reference group contains at least 10 spots to reduce potential bias.
#'
#' @param object TUSCAN object.
#' @param seed Optional numerical value. Set a random seed for reproducible results.
#' @param method Clustering method. Options: c('louvain',"km")
#' @param res Optional resolution parameter in louvain clustering that allows the user to adjust the number of clusters. Larger values typically yield more clusters.
#' @param k Number of clusters from K-Means.
#'
#' @return Return a processed TUSCAN object adding preliminary clustering results.
#' @import dplyr
#' @importFrom FNN get.knn
#' @importFrom igraph graph_from_data_frame simplify cluster_louvain membership
#' @export
#'
#' @examples

PreCluster <- function(object,res=NULL,method=c("louvain","km"),k=5,seed=NULL){

  method <- match.arg(method)

  if(is.null(seed)){
    seed=as.integer(Sys.time())
  }

  img <- object@image
  encoded_data <- object@autoencoder

  global_min <- mean(apply(encoded_data, MARGIN = 2, FUN = min))
  global_max <- mean(apply(encoded_data, MARGIN = 2, FUN = max))

  img <- img %>% mutate(RS_r=rescale_to_range(img$RS,global_min,global_max),
                        GS_r=rescale_to_range(img$GS,global_min,global_max),
                        BS_r=rescale_to_range(img$BS,global_min,global_max))
  features <- cbind(encoded_data,img[,c("RS_r","GS_r","BS_r")])

  if (method == "louvain") {

    set.seed(seed)
    info.spatial = as.data.frame(features)
    colnames(info.spatial) =  paste0("factor", 1:ncol(features))
    knn.norm = get.knn(as.matrix(features), k = 100)
    knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),
                                     k=100), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
    nw.norm = graph_from_data_frame(knn.norm, directed = FALSE)
    nw.norm = simplify(nw.norm)

    if (is.null(res)) {
      lc.norm = cluster_louvain(nw.norm)
    } else {
      lc.norm = cluster_louvain(nw.norm, resolution = res)
    }
    cluster_ids <- lc.norm$membership

  }

  if (method == "km") {

    set.seed(seed)
    cluster_result = kmeans(as.matrix(features), centers = k)
    cluster_ids <- cluster_result$cluster
  }

  print(paste("Found", length(unique(cluster_ids)), "clusters"))
  #print(paste("Found", length(unique(membership(lc.norm))), "clusters"))
  #img$cluster <- factor(lc.norm$membership)
  #object@pre_cluster <- factor(lc.norm$membership)
  img$cluster <- factor(cluster_ids)
  object@pre_cluster <- factor(cluster_ids)
  object@image <- img

  # save clustering results
  df <- data.frame(barcode=img$barcode, cluster=object@pre_cluster)
  row.names(df) <- df$barcode
  object@annotation <- df
  write.table(df, file = "annotations.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

  return(object)
}
