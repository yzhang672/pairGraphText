
# Learning to use the R package: pairGraphText
## Yilin Zhang

### Introduction

The pairGraphText package provides R functions to implement the spectral contextualization method, pairGraphText, when analyzing networks with node covariates.

### Installation

```{r}
install.packages("devtools")
library(devtools)
install_github("yzhang672/pairGraphText")
library(pairGraphText)
library("irlba"); library("ggplot2")
```

### Your First pairGraphText Adventure

#### Simulate a network with node covariates.
```{r}
set.seed(12017)

source("https://raw.githubusercontent.com/karlrohe/fastRG/master/fastRDPG.R")

n = 1000  # number of nodes
d = 1000 # number of covariates
K = 5  # number of blocks

# Here are the parameters for the graph:

pi = rexp(K) +1
pi = pi/sum(pi) * 3
B = matrix(rexp(K^2)+1, nrow=K)
diag(B) = diag(B)+ mean(B)*K
theta = rexp(n)
paraG = dcsbm(theta=theta, pi = pi, B=B,parametersOnly = T)


# Here are the parameters for the covariates:

thetaY = rexp(d)
piFeatures = rexp(K) +1
piFeatures = piFeatures/sum(piFeatures) * 3
BFeatures = matrix(rexp(K^2)+1, nrow=K)
diag(BFeatures) = diag(BFeatures)+ mean(BFeatures)*K

paraFeat = dcsbm(theta = thetaY,pi = piFeatures, B = BFeatures,parametersOnly = T)

# the node "degrees" in the features, should be related to their degrees in the graph.
X = paraG$X
X@x = paraG$X@x + rexp(n) # the degree parameter should be different. X@x + rexp(n) makes feature degrees correlated to degrees in graph.


# Generate the AdjMat and BowMat
AdjMat = fastRG(paraG$X,paraG$S, avgDeg = 10)
BowMat = fastRG(X,paraFeat$S, paraFeat$X, avgDeg = 20)
BowMatX = BowMatY = BowMat;
```

#### Implement PairGraphText

```{r}
source('~/Dropbox/my project/RPackage/pairGraphText/R/pairGraphText.R')
# Set parameter values
signif_level = 0.05; weight_h = 0;
ncluster_x = 4; ncluster_y = 5; niteration = 1000;

# Do pairGraphText
ptm0 <- proc.time()

ClusterResult <-
  pairGraphText(AdjMat, BowMatX, BowMatY, signif_level,
         weight_h, ncluster_x, ncluster_y, niteration)

ptm1 <- proc.time()
print(ptm1-ptm0)
```

#### Visualize your clustering results.

```{r}
source('~/Dropbox/my project/RPackage/pairGraphText/R/SBM.R')
cluster_x <- ClusterResult[[2]]$cluster
cluster_y <- ClusterResult[[3]]$cluster
SbmMat <- createB(AdjMat, cluster_x, cluster_y)

balloonGgPlot(SbmMat, nscale = 1, logTran = FALSE, sqrtTran = FALSE, xlabel = "X Clusters", ylabel = "Y Clusters", main = "SBM Matrix") + theme(plot.title = element_text(hjust = 0.5))
```
