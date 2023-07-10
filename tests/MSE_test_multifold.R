library(countsplit)

overdisp <- 3
K <- 3
clustTry = 2:10
B1 = 1.5
n=400
p=80
ntrials <- 50
folds=5

poisMSEs <- matrix(NA, nrow=ntrials, ncol=length(clustTry))
knownMSEs <- matrix(NA, nrow=ntrials, ncol=length(clustTry))
sctMSEs <- matrix(NA, nrow=ntrials, ncol=length(clustTry))

cluster.sse.log <- function(trainDat, testDat, clusters.train, clusters.test,
                            eps.train, eps.test) {
  totSS <- 0
  for (lab in unique(clusters.train)) {
    if (sum(clusters.test==lab) > 1 &sum(clusters.train==lab) > 1 ) {
    clustdat.test <- testDat[clusters.test==lab,]
    clustdat.train <- trainDat[clusters.train==lab,]
    #### This is done on the scale of the original data X. 
    ### BC of sparsity
    colmeansTrain <- apply(clustdat.train, 2, mean)
    pred.means <- 1/eps.train*  colmeansTrain 
    ss <- apply(1/eps.test*clustdat.test, 1,  function(u) sum((log(u+1)-log(pred.means+1))^2))
    totSS <- totSS+sum(ss)
    } 
  }
  return(totSS)
}

for (trial in 1:ntrials) {
    print(trial)
    set.seed(trial)
   
    B0s <- rnorm(p)
    clusters <- sample(1:K, size=n, replace=TRUE)
    logLambda <- matrix(B0s, nrow=n, ncol=p, byrow = TRUE)
  
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

    X <- sapply(1:p, function(u) rnbinom(length(Lambda[,u]), mu=Lambda[,u], size= true.overdisps[u]))
    partition.pois <- multisplit(X, folds=5)
    partition.known <- multisplit(X, folds=5, overdisps=true.overdisps)
    
    rownames(X) <- 1:NROW(X)
    colnames(X) <- 1:NCOL(X)
    fit <- sctransform::vst(t(X), verbosity=0, min_cells=0)
    overdisps.hat <- fit$model_pars_fit[,1]
    if (length(overdisps.hat) != NCOL(X)) {
      temp <- overdisps.hat
      overdisps.hat<- rep(Inf, NCOL(X))
      overdisps.hat[rownames(temp)] <- temp
    }
    
    partition.hat <- multisplit(X, folds=5, overdisps=overdisps.hat)
    
    ### POISSON
    fullRes <- matrix(NA, nrow=folds, ncol=length(clustTry))
    for (fold in 1:folds) {
        testDat <- partition.pois[[fold]]
        trainDat <- X - testDat
        for (j in 1:length(clustTry)) {
            clusters.train <- kmeans(log(trainDat+1), centers=clustTry[j], nstart=10)$cluster
            fullRes[fold,j] <- cluster.sse.log(trainDat, testDat, clusters.train, clusters.train, 1-1/folds, 1/folds)
        }
    }
    poisMSEs[trial, ] <- colMeans(fullRes)
    
    ### KNOWN
    fullRes <- matrix(NA, nrow=folds, ncol=length(clustTry))
    for (fold in 1:folds) {
      testDat <- partition.known[[fold]]
      trainDat <- X - testDat
      for (j in 1:length(clustTry)) {
        clusters.train <- kmeans(log(trainDat+1), centers=clustTry[j], nstart=10)$cluster
        fullRes[fold,j] <- cluster.sse.log(trainDat, testDat, clusters.train, clusters.train, 1-1/folds, 1/folds)
      }
    }
    knownMSEs[trial, ] <- colMeans(fullRes)
    
    
    ### EST
    fullRes <- matrix(NA, nrow=folds, ncol=length(clustTry))
    for (fold in 1:folds) {
      testDat <- partition.hat[[fold]]
      trainDat <- X - testDat
      for (j in 1:length(clustTry)) {
        clusters.train <- kmeans(log(trainDat+1), centers=clustTry[j], nstart=10)$cluster
        fullRes[fold,j] <- cluster.sse.log(trainDat, testDat, clusters.train, clusters.train, 1-1/folds, 1/folds)
      }
    }
    sctMSEs[trial, ] <- colMeans(fullRes)
}
        
  

logMeanPois <- log(colMeans(poisMSEs))
logMeanKnown <- log(colMeans(knownMSEs))
logMeanSCT <- log(colMeans(sctMSEs))

plot(clustTry, (logMeanPois-min(logMeanPois))/(max(logMeanPois)-min(logMeanPois)), type='b', col="red")
points(clustTry,  (logMeanKnown-min(logMeanKnown))/(max(logMeanKnown)-min(logMeanKnown)), type='b', col="green")
points(clustTry,  (logMeanSCT-min(logMeanSCT))/(max(logMeanSCT)-min(logMeanSCT)), type='b', col="blue")
