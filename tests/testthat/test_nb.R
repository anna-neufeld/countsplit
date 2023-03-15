# library(ggplot2)
# library(countsplit)
# 
# 
# ### Parameters for simulation. 
# n <- 1000
# p <- 1000
# B1 <- 2
# overdisp <- 5
# K <- 3
# 
# ### Generate a dataset!
# set.seed(1)
# B0s <- rnorm(p)
# clusters <- sample(1:K, size=n, replace=TRUE)
# logLambda <- matrix(B0s, nrow=n, ncol=p, byrow = TRUE)
# c <- 1
# if (K > 1) {
#   for (clust in 1:(K-1)) {
#     logLambda[clusters==clust,c:(c+p/20-1)] <-  logLambda[clusters==clust,c:(c+p/20-1)]+B1
#     c <- c+p/20
#   }
# }
# Lambda <- exp(logLambda)
# Lambda_bar_js <- colMeans(Lambda)
# overdisps <- sapply(Lambda_bar_js, function(u) u/overdisp)
# X <- sapply(1:p, function(u) rnbinom(length(Lambda[,u]), mu=Lambda[,u], size=overdisps[u]))
#   
# 
# #### Do count splitting to estimate the number of clusters! Here we assume that the
# #### overdispersion values are unknown and so we use sctransform. 
# 
# ### First, I wrote my own sctransform wrapper function. 
# get_sctransform_bs <- function(X) {
#   rownames(X) <- 1:NROW(X)
#   colnames(X) <- 1:NCOL(X)
#   
#   #### This uses all of sctransform's default settings. 
#   fit <- sctransform::vst(t(X), verbosity=0, min_cells=0)
#   
#   #### I want to make sure that it returns me an estimated overdispersion for EVERY gene.
#   #### by default, it drops the ones that are estimated to be Poisson.
#   #### I add them back in, with an overdisp of infinity.
#   b_j1 <- fit$model_pars_fit[,1]
#   if (length(b_j1) != NCOL(X)) {
#     temp <- b_j1
#     b_j1 <- rep(Inf, NCOL(X))
#     b_j1[rownames(temp)] <- temp
#   }
#   return(b_j1)
# }
# 
# overdisps.sct <- get_sctransform_bs(X)
# 
# ### Sanity check: are these estimates good?
# nulls <- c(rep(FALSE, K*p/20), rep(TRUE, p-K*p/20))
# plot(overdisps.sct, overdisps, col=as.factor(nulls))
# abline(0,1,col="red")
# 
# ### Do the splitting
# eps <- 0.5
# Xsplit <- nb.countsplit(X, epsilon=eps, overdisps=overdisps.sct)
# Xtrain <- Xsplit$train
# Xtest <- Xsplit$test
# 
# #### Sanity check-- are our train/test sets mostly independent?
# cors <- sapply(1:p, function(u) cor(Xtrain[,u], Xtest[,u]))
# plot(cors)
# ### For the null genes, 75% have correlation less than 0.02. That's good! Yay. 
# summary(abs(cors[(K*p/20+1):p]))
# ### For the non-null genes, the correlations are a little bigger but not THAT much bigger. Huh.
# ### Weak signal strength??
# summary(abs(cors[1:(K*p/20)]))
# 
# 
# 
# ### Now let's look at if count splitting can estimate the number of clusters. 
# cluster.sse.log <- function(trainDat, testDat, clusters, eps.train) {
#   eps.test <- 1-eps.train
#   totSS <- 0
#   for (lab in unique(clusters.train)) {
#     clustdat.test <- testDat[clusters==lab,, drop='F']
#     clustdat.train <- trainDat[clusters==lab,, drop='F']
#     
#     #### This is done on the scale of the original data X. 
#     pred.means <- 1/eps.train*colMeans(clustdat.train)
#     ss <- apply(1/eps.test*clustdat.test, 1,  function(u) sum((log(u+1)-log(pred.means+1))^2))
#     totSS <- totSS+sum(ss)
#   }
#   return(totSS)
# }
# 
# clustTry = 10
# rands <- rep(NA, clustTry)
# res <- rep(NA, clustTry)
# for (j in 1:clustTry) {
#   clusters.train <- kmeans(log(Xtrain+1), centers=j, nstart=10)$cluster
#   rands[j] <- mclust::adjustedRandIndex(clusters.train, clusters)
#   res[j] <- cluster.sse.log(Xtrain, Xtest, clusters.train, eps)
# }
# 
# ggplot(data=NULL)+
#   geom_line(aes(x=1:clustTry, y=res))+
#   geom_vline(xintercept=K, col="red")
# 
# #### It helped that the signal was strong enough here to identify almost exactly correct clusters when K=3. 
# rands[K]
# 
