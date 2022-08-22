## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TDAvec)
library(TDA) # to compute persistence diagrams
library(microbenchmark) # to compare computational costs

## -----------------------------------------------------------------------------
N <- 100 # point cloud size
set.seed(123)
X <- circleUnif(N) + rnorm(2*N,mean = 0,sd = 0.2)
# plot the point cloud
plot(X,pch = 20,asp = 1)

## -----------------------------------------------------------------------------
D <- ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram
sum(D[,1]==0) # number of connected components
sum(D[,1]==1) # number of loops
sum(D[,1]==2) # number of voids

## -----------------------------------------------------------------------------
plot(D)
# the solid dots represent connected components
# the red triangles represent loops

## -----------------------------------------------------------------------------
# sequence of scale values to vectorize the summary function
scaleSeq = seq(0,2,length.out=11) 
# compute the PL and PS summaries for homological dimension H_0
computePL(D,homDim = 0,scaleSeq,k=1)
computePS(D,homDim = 0,scaleSeq,p=1)

## -----------------------------------------------------------------------------
# compute the PL and PS summaries for homological dimension H_1
computePL(D,homDim = 1,scaleSeq,k=1)
computePS(D,homDim = 1,scaleSeq,p=1)

## -----------------------------------------------------------------------------
pl1 <- computePL(D,homDim = 0,k=1,scaleSeq)
pl2 <- as.vector(landscape(D,dimension = 0,KK = 1, tseq = scaleSeq))
all.equal(pl1,pl2) # -> TRUE (the results are the same)

compCost <- microbenchmark(
  computePL(D,homDim = 0,k=1,scaleSeq),
  landscape(D,dimension = 0,KK = 1, tseq = scaleSeq),
  times = 500
)
sm <- summary(compCost)
costRatioPL <- sm$mean[2]/sm$mean[1] # ratio of computational time means
print(costRatioPL)

## -----------------------------------------------------------------------------
ps1 <- computePS(D,homDim = 0, p = 1,scaleSeq)
ps2 <- as.vector(silhouette(D,dimension = 0,p = 1, tseq = scaleSeq))
all.equal(ps1,ps2) # -> TRUE (the results are the same)

compCost <- microbenchmark(
  computePS(D,homDim = 0, p = 1,scaleSeq),
  silhouette(D,dimension = 0,p = 1, tseq = scaleSeq),
  times = 500
)
sm <- summary(compCost)
costRatioPS <- sm$mean[2]/sm$mean[1]
print(costRatioPS)

## -----------------------------------------------------------------------------
# Persistent Entropy Summary (PES) function
# compute PES for homological dimension H0
computePES(D,homDim = 0,scaleSeq) 
# compute PES for homological dimension H1
computePES(D,homDim = 1,scaleSeq)

# Euler Characteristic Curve (ECC) 
computeECC(D,maxhomDim = 1,scaleSeq) # maxhomDim = maximal homological dimension considered

# Vector of Averaged Bettis (VAB) - a vectorization of Betti Curve
# compute VAB for homological dimension H0
computeVAB(D,homDim = 0,scaleSeq) 
# compute VAB for homological dimension H1
computeVAB(D,homDim = 1,scaleSeq)

## -----------------------------------------------------------------------------
D[,3] <- D[,3] - D[,2] 
colnames(D)[3] <- "Persistence"

## -----------------------------------------------------------------------------
# Persistence Image (PI)
res <- 5 # resolution or grid size
# find min and max persistence values
minPH0 <- min(D[D[,1]==0,3]); maxPH0 <- max(D[D[,1]==0,3])
sigma <- 0.5*(maxPH0-minPH0)/res # default way of selecting the standard deviation sigma of the Gaussians on top of each point of the diagram
# construct one-dimensional grid of scale values
ySeqH0 <- seq(minPH0,maxPH0,length.out=res+1)
# compute PI for homological dimension H_0
computePI(D,homDim=0,xSeq=NA,ySeqH0,res,sigma)

# Vectorized Persistence Block (VPB)
# construct one-dimensional grid of scale values
ySeqH0 <- unique(quantile(D[D[,1]==0,3],probs = seq(0,1,by=0.2))) 
tau <- 0.3 # parameter in [0,1] which controls the size of blocks around each point of the diagram 
# compute VPB for homological dimension H_0
computeVPB(D,homDim = 0,xSeq=NA,ySeqH0,tau) 

## -----------------------------------------------------------------------------
# PI
res <- 5 # resolution or grid size
# find min & max birth and persistence values
minBH1 <- min(D[D[,1]==1,2]); maxBH1 <- max(D[D[,1]==1,2])
minPH1 <- min(D[D[,1]==1,3]); maxPH1 <- max(D[D[,1]==1,3])
xSeqH1 <- seq(minBH1,maxBH1,length.out=res+1)
ySeqH1 <- seq(minPH1,maxPH1,length.out=res+1)
sigma <- 0.5*(maxPH1-minPH1)/res
# compute PI for homological dimension H_1
computePI(D,homDim=1,xSeqH1,ySeqH1,res,sigma)

# VPB
xSeqH1 <- unique(quantile(D[D[,1]==1,2],probs = seq(0,1,by=0.2)))
ySeqH1 <- unique(quantile(D[D[,1]==1,3],probs = seq(0,1,by=0.2)))
tau <- 0.3
# compute VPB for homological dimension H_1
computeVPB(D,homDim = 1,xSeqH1,ySeqH1,tau) 

