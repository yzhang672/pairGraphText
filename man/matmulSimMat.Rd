% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpeCon.R
\name{matmulSimMat}
\alias{matmulSimMat}
\title{matmulSimMat
Self-defined matrix multiplication function for fast SVD on the Similarity Matrix}
\usage{
matmulSimMat(GraphLap, weight_h, Xs, Ys, CallRespMat_thresh)
}
\arguments{
\item{GraphLap}{regularized Graph Laplacian}

\item{weight_h}{weight of the text assisted part h}

\item{Xs}{Scaled Bag of Word Matrix for the first type of nodes (senders)}

\item{Ys}{Scaled Bag of Word Matrix for the second type of nodes (receivers)}

\item{CallRespMat_thresh}{Call Response Matrix after Thresholding}
}
\description{
matmulSimMat
Self-defined matrix multiplication function for fast SVD on the Similarity Matrix
}
\keyword{SVD}
\keyword{irlba}
\keyword{matmul}
