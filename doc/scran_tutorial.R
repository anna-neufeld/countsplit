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

## -----------------------------------------------------------------------------
cm

## -----------------------------------------------------------------------------
dim(counts(cm))
counts(cm)[1:10,1:10]

## -----------------------------------------------------------------------------
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

