
#' Predict tumor spots
#'
#' This function implements the automatic tumor segmentation and classification.
#' @param object TUSCAN object.
#' @param seed Optional numerical value. Set a random seed for reproducible results.
#' @param maxK Maximum cluster number to evaluate.
#' @param reps Number of subsamples.
#' @param pItem Proportion of items to sample.
#' @param pFeature Proportion of features to sample.
#' @param title File directory for storing consensus clustering figures.
#' @param clusterAlg Cluster algorithm. (default: 'km'), options('km','hc','pam'). 'km': k-means, 'hc': hierarchical (hclust), 'pam': paritioning around medoids.
#' @param distance Distance function. (default: 'euclidean'), options('euclidean','pearson','spearman','maximum','canberra','minkowski').
#'
#' @return Return a a processed TUSCAN object adding tumor annotations.
#' @import ConsensusClusterPlus
#' @export
#'
#' @examples

TumorCluster <- function(object,seed=NULL,maxK=3,reps=100,pItem=0.8,pFeature=1,
                         title=NULL,clusterAlg=c("km","hc","pam"),
                         distance=c("euclidean","pearson","spearman","maximum","canberra","minkowski")){

  clusterAlg <- match.arg(clusterAlg)
  distance <- match.arg(distance)

  if(is.null(title)){
    title <- tempdir()
  }

  if(is.null(seed)){
    seed=as.integer(Sys.time())
  }

  d <- object@cnv.data
  results = ConsensusClusterPlus(d,maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,
                                 title=title,clusterAlg=clusterAlg,distance=distance,seed=seed,plot="png")

  #classify all spots into two clusters by consensus clustering
  consensus_label <- results[[2]][["consensusClass"]]

  img <- object@image
  img <- img[match(names(consensus_label),img$barcode),]
  img$consensus <- factor(consensus_label)

  normal_loc <- unlist(object@reference_grouped_spot_indices)
  norm.mat <- t(d[,normal_loc])
  id1 <- img %>% filter(consensus==1) %>% select(barcode) %>% pull(1)
  id2 <- img %>% filter(consensus==2) %>% select(barcode) %>% pull(1)
  mat1 <- t(d[,id1])
  mat2 <- t(d[,id2])

  ##assign labels to two clusters
  #calculate the distance between the centroids of two clusters and that of the normal reference.
  #the cluster that is further away from the norm is annotated as tumor.
  centroids <- lapply(list(norm.mat, mat1, mat2), function(cluster) {
    colMeans(cluster, na.rm = TRUE)
  })
  centroids_matrix <- do.call(rbind, centroids)
  rownames(centroids_matrix) <- c("Norm", "Mat1", "Mat2")
  dist_matrix <- proxy::dist(centroids_matrix)
  tumor <- which.max(dist_matrix[1:2])
  img <- img %>% mutate(tumor = ifelse(consensus == tumor, "Tumor", "Normal"),
                        tumor = factor(tumor, levels = c("Tumor", "Normal")))
  img <- img[match(colnames(object@cnv.data),img$barcode),]
  object@tumor_cluster <- img$tumor

  object@image <- img

  return(object)

}
