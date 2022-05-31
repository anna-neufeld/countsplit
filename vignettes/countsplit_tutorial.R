## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE
)

## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github("anna-neufeld/countsplit")

## ---- message=FALSE-----------------------------------------------------------
library(countsplit)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 1000
p <- 200
X <- matrix(rpois(n*p, lambda=5), nrow=n)

## -----------------------------------------------------------------------------
set.seed(3)
clusters.full <- kmeans(log(X+1), centers=2)$cluster
results.naive <- fit_glms(X, clusters.full, family="poisson")
head(results.naive)

## -----------------------------------------------------------------------------
library(ggplot2)
ggplot(data=results.naive, aes(sample=pval))+geom_qq(distribution="qunif")+geom_abline(col="red")

## -----------------------------------------------------------------------------
set.seed(2)
split <- countsplit(X, epsilon=0.5)
Xtrain <- split$train
Xtest <- split$test
clusters.train <- kmeans(log(Xtrain+1), centers=2)$cluster
results.countsplit <- fit_glms(Xtest, clusters.train, family="poisson")
ggplot(data=results.countsplit,aes(sample=pval))+geom_qq(distribution="qunif")+geom_abline(col="red")

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(scran)
library(tidyverse)

## -----------------------------------------------------------------------------
data(cm)

## -----------------------------------------------------------------------------
set.seed(1)
X <- t(counts(cm))
split <- countsplit(X, epsilon=0.5)
Xtrain <- split$train
Xtest <- split$test

## -----------------------------------------------------------------------------
cm.train <- cm
counts(cm.train) <- t(Xtrain)

## -----------------------------------------------------------------------------
sizeFactors(cm.train) <- librarySizeFactors(cm.train)
## Pre-processing
cm.train <- logNormCounts(cm.train)
top.hvgs <- getTopHVGs(modelGeneVar(cm.train), fdr.threshold=0.05)
cm.train <- fixedPCA(cm.train, subset.row=top.hvgs)

## Graph clustering
clusters.train <- clusterCells(cm.train,use.dimred="PCA")

## ----out.width="90%"----------------------------------------------------------
table(clusters.train)
ggplot(as_tibble(reducedDim(cm.train)), aes(x=PC1, y=PC2, col=as.factor(clusters.train)))+geom_point()+labs(col="Cluster")

## ---- warning=FALSE-----------------------------------------------------------
set.seed(1)
indices <- which(clusters.train==1 | clusters.train==2)
genes <- sample(1:NCOL(Xtest), size=500)
results <- fit_glms(Xtest[indices,genes], clusters.train[indices], offsets=clusters.train$sizeFactors, family="quasipoisson")
table(results$pval < 0.01)
head(results)

## ---- warning=FALSE-----------------------------------------------------------
library(Seurat)

## ---- eval=FALSE--------------------------------------------------------------
#  library(scran)
#  data(cm)
#  sizeFactors(cm) <- librarySizeFactors(cm)
#  cm <- logNormCounts(cm)
#  cm.seurat <- as.Seurat(cm)

