#' createB
#'
#' Create SBM matrix
#' @param AdjMat An Adjacency Matrix
#' @param kms_cluster A vector of cluster assignments for the first type of nodes (senders)
#' @param kmr_cluster A vector of cluster assignments for the second type of nodes (receivers)
#'
#' @return SbmMat is the SBM matrix
#' @keywords SBM Matrix
#' @export

createB <- function(AdjMat, kms_cluster, kmr_cluster)
{
  Zs_hat = model.matrix(~as.factor(kms_cluster)-1);
  Zr_hat = model.matrix(~as.factor(kmr_cluster)-1);
  zzs = solve(t(Zs_hat)%*%Zs_hat); zzr = solve(t(Zr_hat)%*%Zr_hat);
  SbmMat = zzs%*%t(Zs_hat)%*%AdjMat%*%Zr_hat%*%zzr
  rownames(SbmMat) = paste("Cluster_x", 1:length(unique(kms_cluster)))
  colnames(SbmMat) = paste("Cluster_y", 1:length(unique(kmr_cluster)))

  return(SbmMat)
}

#' balloonGgPlot
#'
#' @param M A matrix
#' @param nscale scaling number
#' @param logTran Whether to do log-transformation
#' @param sqrtTran Whether to do sqrt-transformation
#' @param xlabel label for rows
#' @param ylabel label for columns
#' @param main title
#'
#' @return A balloon plot for the matrix
#' @keywords BalloonPlot ggplot
#'
#' @export

balloonGgPlot = function(M, nscale, logTran, sqrtTran, xlabel, ylabel, main){

  M = t(M)
  n = nrow(M)
  d = ncol(M)

  M = M[,d:1]*nscale;

  scaleItMean = mean(abs(M))

  if(logTran){
    scaleItMean = mean(log(abs(M) +1))
  }
  if(sqrtTran){
    scaleItMean = mean(sqrt(abs(M)))
  }

  M = as.data.frame(as.table(as.matrix(M)))

  DotSize = rep(0,n*d);

  for(i in 1:(n*d)){
    dotSize = abs(M$Freq[i])
    if(logTran) dotSize= log(dotSize + 1)
    if(sqrtTran) dotSize = sqrt(dotSize)
    DotSize[i] = dotSize

  }

  M = data.frame(M, DotSize)
  colnames(M)[c(1,2)] = c("Cluster_y_ID", "Cluster_x_ID")

  s <- ggplot(M, aes(Cluster_y_ID, Cluster_x_ID)) +
    geom_point( aes(size = DotSize))+
    xlab(ylabel) +
    ylab(xlabel) +
    ggtitle(main) +
    theme_bw() +
    theme(
      #legend.position = "none",
      axis.text.x = element_text(size=10, face = "bold"),
      axis.text.y = element_text(size=10, face = "bold"),
      axis.title.x = element_text(size=12, face = "bold"),
      axis.title.y = element_text(size=12, face = "bold"),
      title = element_text(size = 14, face = "bold")
    )
  s
}
