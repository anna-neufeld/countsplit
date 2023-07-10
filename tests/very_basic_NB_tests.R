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
true.overdisps <- sapply(Lambda_bar_js, function(u) u/overdisp)
X <- sapply(1:p, function(u) rnbinom(length(Lambda[,u]), mu=Lambda[,u], size=true.overdisps[u]))


# Try Poisson count splitting.
poissonSplit <- nb.countsplit(X)
xTrain_pois <- poissonSplit$train
xTest_pois <- poissonSplit$test
cors <- sapply(1:p, function(u) cor(xTrain_pois[,u], xTest_pois[,u]))
mean(cors[1:10])
mean(cors[11:p])
mclust::adjustedRandIndex(kmeans(log(xTrain_pois+1), 2)$cluster, kmeans(log(xTest_pois+1), 2)$cluster)

# Try negative binomial count splitting with known overdispersion. 
nbcsSplit_known <- nb.countsplit(X, overdisps = overdisps)
xTrain_known <- nbcsSplit_known$train
xTest_known <- nbcsSplit_known$test
cors <- sapply(1:p, function(u) cor(xTrain_known[,u], xTest_known[,u]))
mean(cors[1:10])
mean(cors[11:p])
mclust::adjustedRandIndex(kmeans(log(xTrain_known+1), 2)$cluster, kmeans(log(xTest_known+1), 2)$cluster)


# Try negative 

##### Sparse matrices