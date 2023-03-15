## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("scran")

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(scran)
library(tidyverse)
library(countsplit)

## -----------------------------------------------------------------------------
data(cm)
cm

## ---- collapse=TRUE-----------------------------------------------------------
dim(counts(cm))
counts(cm)[1:10,1:10]

## ---- collapse=TRUE-----------------------------------------------------------
table(cm$individual)
table(cm$diffday)

## -----------------------------------------------------------------------------
set.seed(1)
X <- counts(cm)
split <- countsplit(X, epsilon=0.5)
Xtrain <- split$train
Xtest <- split$test

## -----------------------------------------------------------------------------
cm.train <- cm
counts(cm.train) <- Xtrain

## -----------------------------------------------------------------------------
clusters <- quickCluster(cm.train)
cm.train <- computeSumFactors(cm.train, clusters=clusters)
cm.train <- logNormCounts(cm.train)

## -----------------------------------------------------------------------------
top.hvgs <- getTopHVGs(modelGeneVar(cm.train), n=2000)

## -----------------------------------------------------------------------------
cm.train <- fixedPCA(cm.train, subset.row=top.hvgs)

## -----------------------------------------------------------------------------
clusters.train <- clusterCells(cm.train,use.dimred="PCA")

## ----out.width="90%", results='hide', warning=FALSE---------------------------
table(clusters.train)
ggplot(as_tibble(reducedDim(cm.train)), aes(x=PC1, y=PC2, col=as.factor(clusters.train)))+geom_point()+labs(col="Cluster")

## ---- warning=FALSE, collapse=TRUE--------------------------------------------
set.seed(1)
indices <- which(clusters.train==1 | clusters.train==2)
genes <- sample(1:NCOL(Xtest), size=500) 
Xtestsubset <- Xtest[genes, indices]
sizefactors.subset <- sizeFactors(cm.train)[indices]

## ---- warning=FALSE, collapse=TRUE--------------------------------------------
results <- t(apply(Xtestsubset, 1, function(u) summary(glm(u~clusters.train[indices], offset=log(sizefactors.subset), family="poisson"))$coefficients[2,]))
table(results[,4] < 0.01) 

## ---- warning=FALSE, collapse=TRUE--------------------------------------------
head(results[results[,4] < 0.01,], n=10)

## -----------------------------------------------------------------------------
cm.test <- SingleCellExperiment(list(counts=Xtest))

## -----------------------------------------------------------------------------
sizeFactors(cm.test) <- sizeFactors(cm.train)
cm.test <- logNormCounts(cm.test) 

## -----------------------------------------------------------------------------
results <- scran::findMarkers( 
cm.test, groups= clusters.train,
 pval.type = "all") 
results[[1]]

