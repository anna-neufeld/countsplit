## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(monocle3)
library(tidyverse)
library(countsplit)

## -----------------------------------------------------------------------------
expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"))

## -----------------------------------------------------------------------------
epsilon=0.5
set.seed(1)
Xtrain <-apply(expression_matrix,2,function(u) rbinom(n=length(u), size=u, p=epsilon))
Xtest <- expression_matrix-Xtrain
cds_train <- new_cell_data_set(Xtrain,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

