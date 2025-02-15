## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TDAvec)
library(TDAstats) # to compute persistence diagrams

## -----------------------------------------------------------------------------
# the number of points to sample
N <- 100 
# set a random seed for reproducibility
set.seed(123) 
# sample N points uniformly from the unit circle and add Gaussian noise
theta <- runif(N, min = 0, max = 2 * pi)
X <- cbind(cos(theta), sin(theta)) + rnorm(2 * N, mean = 0, sd = 0.2)
# plot the point cloud
plot(X,pch = 20,asp = 1,xlab = 'x',ylab = 'y')

## -----------------------------------------------------------------------------
D <- calculate_homology(X,dim=1,threshold=2)
sum(D[,1]==0) # number of connected components
sum(D[,1]==1) # number of loops
sum(D[,1]==2) # number of voids

## -----------------------------------------------------------------------------
plot_persist(D)

## -----------------------------------------------------------------------------
# sequence of scale values to be used for vectorization
scaleSeq = seq(0,1.5,length.out=16) 
# vectorize the Betti curve for homological dimension H_1
computeBettiCurve(D,homDim=1,scaleSeq)

## -----------------------------------------------------------------------------
# classic
computePersistenceLandscape(D,homDim=1,scaleSeq,k=3) # k = the number of persistence landscape functions to consider (default is 1)

# generalized
computePersistenceLandscape(D,homDim=1,scaleSeq,k=3,generalized=TRUE,kernel = "epanechnikov",h=0.2) # h = bandwidth for the kernel function

## -----------------------------------------------------------------------------
computeAlgebraicFunctions(D,homDim=0)  

## -----------------------------------------------------------------------------
D[,3] <- D[,3] - D[,2] 
colnames(D)[3] <- "Persistence"

## -----------------------------------------------------------------------------
# Persistence Image (PI)
resB <- 5 # resolution (or grid size) along the birth axis
resP <- 5 # resolution (or grid size) along the persistence axis 
# find min and max persistence values for dimension H_0
minPH0 <- min(D[D[,1]==0,3]); maxPH0 <- max(D[D[,1]==0,3])
# construct one-dimensional grid of scale values
ySeqH0 <- seq(minPH0,maxPH0,length.out=resP+1)
# default way of selecting the standard deviation sigma of the Gaussians on top of each point of the diagram
sigma <- 0.5*(maxPH0-minPH0)/resP 
# compute PI for homological dimension H_0
computePersistenceImage(D,homDim=0,xSeq=NA,ySeqH0,sigma)

## -----------------------------------------------------------------------------
# Persistence Image (PI)
# find min & max birth and persistence values for dimension H_1
minBH1 <- min(D[D[,1]==1,2]); maxBH1 <- max(D[D[,1]==1,2])
minPH1 <- min(D[D[,1]==1,3]); maxPH1 <- max(D[D[,1]==1,3])
xSeqH1 <- seq(minBH1,maxBH1,length.out=resB+1)
ySeqH1 <- seq(minPH1,maxPH1,length.out=resP+1)
sigma <- 0.5*(maxPH1-minPH1)/resP
# compute PI for homological dimension H_1
computePersistenceImage(D,homDim=1,xSeqH1,ySeqH1,sigma) 

