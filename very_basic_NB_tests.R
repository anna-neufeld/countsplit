##### Regular matrices

# Generate data with two true clusters. 
n=1000
p=1000
B1 = 2
K = 2
overdisp = 5
B0s <- rnorm(p)
clusters <- sample(1:K, size=n, replace=TRUE)
logLambda <- matrix(B0s, nrow=n, ncol=p, byrow = TRUE)
K <- length(unique(clusters))
c <- 1
if (K > 1) {
  for (clust in 1:(K-1)) {
    logLambda[clusters==clust,c:(c+p/20-1)] <-  logLambda[clusters==clust,c:(c+p/20-1)]+B1
    c <- c+p/20
  }
}
Lambda <- exp(logLambda)
Lambda_bar_js <- colMeans(Lambda)
overdisps <- sapply(Lambda_bar_js, function(u) constant_b_function(u, overdisp))
X <- sapply(1:p, function(u) rnbinom(length(Lambda[,u]), mu=Lambda[,u], size=overdisps[u]))


# Try Poisson count splitting, NBCS known, and NBCS estimated

poissonSplit <- countsplit(X)





##### Sparse matrices