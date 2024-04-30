#' Find a subset of normal spots
#'
#' This function helps determine the cluster with the highest potential purity of normal spots through an evaluation score \eqn{s} that combines both gene expression and histology image information.
#' The evaluation score \eqn{s} of each cluster is calculated as \eqn{s=w \times g - (1-w) \times \sigma}, where \eqn{g} is the median gray value of all spots within that cluster extracted from the gray-scaled histology image, \eqn{\sigma} is standard deviation of all gene expression levels within that cluster, and \eqn{w} is the weight of image. Both \eqn{g} and \eqn{\sigma} are standardized.
#' The cluster with the largest \eqn{s} is selected as the normal reference.
#'
#' @param object TUSCAN object.
#' @param w Weight of image. Default value is 0.5, i.e., gene expression and image are equally weighted.
#'
#' @return Print the name of the selected normal cluster.
#' @export
#'
#' @examples
#'

FindNormalCluster <- function(object,w=0.5){

  count.mat <- object@counts.data
  norm.mat <- log(sqrt(count.mat)+sqrt(count.mat+1))
  norm.mat.center <- apply(norm.mat,2,function(x)(x <- x-mean(x)))
  colnames(norm.mat.center) <-  colnames(count.mat)

  ##smooth data
  dlm.sm <- function(c){
    model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
    x <- dlm::dlmSmooth(norm.mat.center[, c], model)$s
    x<- x[2:length(x)]
    x <- x-mean(x)
  }

  test.mc <-parallel::mclapply(1:ncol(norm.mat.center), dlm.sm)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat.center), byrow = FALSE)
  colnames(norm.mat.smooth) <- colnames(norm.mat.center)

  clst <- object@pre_cluster
  group <- unique(clst)

  SDM <-NULL
  SSD <-NULL

  ##estimate gene expression variances across all clusters via GMM model
  for(i in 1:length(group)){
    data.c <- apply(norm.mat.smooth[, which(clst==group[i])],1, median)
    sx <- max(c(0.05, 0.5*sd(data.c)))
    GM3 <- mixtools::normalmixEM(data.c, lambda = rep(1,3)/3, mu = c(-0.2, 0, 0.2), sigma = sx,arbvar=FALSE,ECM=FALSE,maxit=5000)
    SDM <- c(SDM, GM3$sigma[1])
    SSD <- c(SSD, sd(data.c))
    i <- i+1
  }

  img <- object@image
  avg_color <- aggregate(gray ~ cluster, data = img, FUN = median)
  g <- avg_color$gray
  sigma <- SDM
  g.s <- scale(g)
  sigma.s <- scale(sigma)
  ##calculate the evaluation score:s
  s <- w*g.s-(1-w)*sigma.s
  ref_group <- group[which.max(s)]

  print(paste("Cluster ",ref_group," is selected as the normal reference.",sep=""))

}
