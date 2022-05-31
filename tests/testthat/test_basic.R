context("basic testing")
library(countsplit)

test_that("Fake Data QQplot", {
  
  set.seed(1)
  n <- 1000
  p <- 200
  X <- matrix(rpois(n*p, lambda=5), nrow=n)
  set.seed(2)
  split <- countsplit(X, epsilon=0.5)
  Xtrain <- split$train
  Xtest <- split$test
  clusters.train <- kmeans(log(Xtrain+1), centers=2)$cluster
  results.countsplit <- fit_glms(Xtest, clusters.train, family="poisson")
  ggplot(data=results.countsplit, aes(sample=pval))+geom_qq(distribution="qunif")+geom_abline(col="red")

  
  expect_true(length(results.countsplit$pval), n)
  expect_true(min(results.countsplit$pval) >= 0 & max(results.countsplit$pval <=1))
  
  }
)

test_that("Real Data", {
  library(scran)
  library(tidyverse)
  data(cm)
  set.seed(1)
  X <- t(counts(cm))
  split <- countsplit(X, epsilon=0.5)
  Xtrain <- split$train
  Xtest <- split$test
  
  cm.train <- cm
  counts(cm.train) <- t(Xtrain)
  
  sizeFactors(cm.train) <- librarySizeFactors(cm.train)
  cm.train <- logNormCounts(cm.train)
  top.hvgs <- getTopHVGs(modelGeneVar(cm.train), fdr.threshold=0.05)
  cm.train <- fixedPCA(cm.train, subset.row=top.hvgs)
  pt.train <- as_tibble(reducedDim(cm.train))$PC1
  
  results.countsplit <- fit_glms(Xtest, pt.train, family="nb")
  expect_true(length(results.countsplit$pval), n)
  expect_true(min(results.countsplit$pval) >= 0 & max(results.countsplit$pval <=1))

})
