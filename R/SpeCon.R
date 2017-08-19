#' SpeCon
#'
#' @param AdjMat An Adjacency Matrix
#' @param BowMatX A Bag of Word Matrix for the first type of nodes (senders)
#' @param BowMatY A Bag of Word Matrix for the second type of nodes (receivers)
#' The two types of nodes can be the same.
#' @param signif_level Significance level that decide the percentage of
#' non-zero elements of the Call Response Matrix we use. (Thresholding Step)
#' @param weight_h The weight h for the covariate-assisted part
#' @param ncluster_x Number of clusters for the first type of nodes (senders)
#' @param ncluster_y Number of clusters for the second type of nodes (receivers)
#' @param niteration Number of iterations in k-means.
#'
#' @return A list of parameters:
#' svdSimMat Singular values and singular vectors of the similarity matrix.
#' km_x Clustering results for the first type of nodes (senders)
#' km_y Clustering results for the second type of nodes (receivers)
#'
#' @keywords Spectral Contextualization
#'
#' @export
SpeCon <- function(
  AdjMat, BowMatX, BowMatY, signif_level, weight_h, ncluster_x, ncluster_y, niteration)
{
  # Calculate the regularized Graph Laplacian from the Adjacency Matrix
  GraphLap <- laplacian(AdjMat, type = "both", regularizer = TRUE)
  # Scale the Bag of Word Matrices.
  # The centering step is woven in the algorithm for computational efficiency.
  # Woven in two places (1. Calculate Call Response Matrix, 2. matmulSimMat)
  Xs <- center_scale(BowMatX, center = FALSE, scale = TRUE)
  Ys <- center_scale(BowMatY, center = FALSE, scale = TRUE)

  # Calculate the Call Response Matrix (equivalent to using centered BowMat)
  CallRespMat_pre_center <- t(Xs) %*% (GraphLap %*% Ys)
  colMeanX <- as.matrix(colMeans(Xs)); colMeanY <- as.matrix(colMeans(Ys));
  CallRespMat <- CallRespMat_pre_center - sum(GraphLap@x)*colMeanX%*%t(colMeanY)

  # Threshold the Call Respond Matrix
  thresh <- quantile(abs(CallRespMat@x), prob = 1-signif_level)
  CallRespMat_thresh <- threshold(CallRespMat,thresh)

  # No need to calculate the Similarity Matrix
  #SimMat <- GraphLap + weight_h*(Xs%*%CallRespMat_thresh%*%t(Ys) -
  #matrix(t(colMeanX)%*%CallRespMat_thresh%*%colMeanY, nrow = nrow(AdjMat), ncol = ncol(AdjMat)))

  # SVD on SimMat
  # This step can be accelerated by sef-defining matmul
  ncluster <- min(ncluster_x, ncluster_y)
  svdSimMat <- irlba(GraphLap, nu = ncluster, nv = ncluster, matmul = matmulSimMat(GraphLap,weight_h,Xs,Ys,CallRespMat_thresh))

  # Normalize singular vectors to have rowsum 1
  u <- svdSimMat$u[,1:ncluster]
  u <- t(apply(u,1,function(x) return(x/sqrt(sum(x^2)))))
  v <- svdSimMat$v[,1:ncluster]
  v <- t(apply(v,1,function(x) return(x/sqrt(sum(x^2)))))

  # K-means on the normalized singular vectors.
  km_x <- kmeans(u, centers = ncluster_x, nstart = niteration)
  km_y <- kmeans(v, centers = ncluster_y, nstart = niteration)

  # km_x and km_y are the clustering results for two types of nodes:
  # senders and receivers
  # (or in other words: nodes on the rows and nodes on the columns)
  return(list(svdSimMat, km_x, km_y))
}


#'----------------------------------------------------------------------
#'Helper functions
#'----------------------------------------------------------------------

#' laplacian
#'
#' The function to calculate Graph Laplacian from an Adjacency Matrix
#'
#' @param AdjMat An Adjacency Matrix
#' @param type Type to normalize to matrix: only by rows, only by columns, or both
#' @param regularizer Whether to regularize: TRUE or FALSE
#'
#' @return GraphLap The Graph Laplacian
#'
#' @keywords Graph Laplacian
#' @export

laplacian <- function(AdjMat,type,regularizer)
{
  rs = rowSums(AdjMat); cs = colSums(AdjMat);
  taur = 0; tauc = 0
  if(regularizer==TRUE)
  {taur = mean(rs); tauc = mean(cs)}

  if(type=="row")
  {
    GraphLap = Diagonal(length(rs), 1/(rs+taur))%*%AdjMat
  }
  if(type=="column")
  {
    GraphLap = AdjMat%*%Diagonal(length(cs), 1/(cs+tauc))
  }
  if(type=="both")
  {
    GraphLap = Diagonal(length(rs), 1/sqrt(rs+taur))%*%
      AdjMat%*%Diagonal(length(cs), 1/sqrt(cs+tauc))
  }
  return(GraphLap)
}

#' center_scale
#'
#' the function to center and scale the Bag of Word Matrix
#' first scale then center, otherwise NaN's can be produced.
#'
#' @param Mat An Matrix
#' @param center Whether to center or not: TRUE or FALSE
#' @param scale Whether to scale or not: TRUE or FALSE
#'
#' @return Mat The Matrix after processing
#'
#' @keywords Center Scale
#'
#' @export


center_scale <- function(Mat, center = FALSE, scale = TRUE)
{
  if(scale){
    Mat <- laplacian(Mat, type = "both", regularizer = FALSE)
  }
  if(center){
    Mat <- scale(Mat, center = TRUE, scale = FALSE)
  }
  return(Mat)
}

#' threshold
#' the function to threshold Matrix
#'
#' @param Mat A matrix
#' @param thresh Threshold
#'
#' @export
#' @return Mat_thresh The matrix after thresholding
#'
#' @keywords Threshold

threshold <- function(Mat,thresh)
{
  Mat_thresh = Mat
  Mat_thresh[abs(Mat) <= thresh] = 0
  return(Mat_thresh)
}

#' matmulSimMat
#' Self-defined matrix multiplication function for fast SVD on the Similarity Matrix
#'
#' @param GraphLap regularized Graph Laplacian
#' @param weight_h weight of the text assisted part h
#' @param Xs Scaled Bag of Word Matrix for the first type of nodes (senders)
#' @param Ys Scaled Bag of Word Matrix for the second type of nodes (receivers)
#' @param CallRespMat_thresh Call Response Matrix after Thresholding
#'
#' @keywords SVD matmul irlba


matmulSimMat <- function(GraphLap,weight_h,Xs,Ys,CallRespMat_thresh)
{
  v <- apply(GraphLap,2,mean)
  function(GraphLap,x,transpose=FALSE)
  {

    if(transpose)
      return( as.matrix(t(crossprod(x,GraphLap)) - sum(x) * v))

    if(is.infinite(weight_h)){
      as.matrix(  Xs %*% (CallRespMat_thresh %*% crossprod(Ys,x))
                  - matrix(t(as.matrix(colMeans(Xs)))%*%CallRespMat_thresh%*%
                             as.matrix(colMeans(Ys)), nrow = nrow(GraphLap), ncol = ncol(GraphLap)) %*% x
                  - cbind(rep(crossprod(v,x)[1],nrow(GraphLap))))
    }else{
      as.matrix(GraphLap %*% x
                + weight_h*Xs %*% (CallRespMat_thresh %*% crossprod(Ys,x))
                - weight_h*matrix(t(as.matrix(colMeans(Xs)))%*%CallRespMat_thresh%*%
                                    as.matrix(colMeans(Ys)), nrow = nrow(GraphLap), ncol = ncol(GraphLap)) %*% x
                - cbind(rep(crossprod(v,x)[1],nrow(GraphLap))))
    }
  }
}



