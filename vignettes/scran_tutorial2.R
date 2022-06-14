## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("scran")

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(scran)
library(tidyverse)
library(countsplit)

## -----------------------------------------------------------------------------
library(scRNAseq)
sce <- GrunPancreasData()
sce

## -----------------------------------------------------------------------------
X <- counts(sce)
set.seed(1)
split <- countsplit(pbmc.counts, epsilon=0.5)
Xtrain <- split$train
Xtest <- split$test

