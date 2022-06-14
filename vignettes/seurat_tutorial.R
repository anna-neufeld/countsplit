## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("Seurat")

## -----------------------------------------------------------------------------
library(Seurat)
library(countsplit)
library(patchwork)
data(pbmc.counts, package="countsplit")

## -----------------------------------------------------------------------------
rownames(pbmc.counts) <- sapply(rownames(pbmc.counts), function(u) stringr::str_replace_all(u, "_","-"))

## -----------------------------------------------------------------------------
set.seed(1)
split <- countsplit(pbmc.counts, epsilon=0.5)
Xtrain <- split$train
Xtest <- split$test

## -----------------------------------------------------------------------------
pbmc.train <- CreateSeuratObject(counts = Xtrain, min.cells = 3, min.features = 200)

## -----------------------------------------------------------------------------
pbmc.train[["percent.mt"]] <- PercentageFeatureSet(pbmc.train, pattern = "^MT-")
pbmc.train <- subset(pbmc.train, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## -----------------------------------------------------------------------------
dim(Xtrain)
dim(Xtest)
dim(pbmc.train)

## -----------------------------------------------------------------------------
rows <- rownames(pbmc.train)
cols <- colnames(pbmc.train)
Xtestsubset <- Xtest[rows,cols]
dim(Xtestsubset)

## -----------------------------------------------------------------------------
all.equal(GetAssayData(pbmc.train, "counts"), GetAssayData(pbmc.train, "data"))

## -----------------------------------------------------------------------------
GetAssayData(pbmc.train, "scale.data")

## -----------------------------------------------------------------------------
pbmc.train <- NormalizeData(pbmc.train)
all.equal(GetAssayData(pbmc.train, "counts"), GetAssayData(pbmc.train, "data"))

## -----------------------------------------------------------------------------
dim(pbmc.train)
pbmc.train <- FindVariableFeatures(pbmc.train, selection.method = "vst", nfeatures = 2000)
dim(pbmc.train)

## ---- out.width="100%"--------------------------------------------------------
library(ggplot2)
top10 <- head(VariableFeatures(pbmc.train), 10)
plot1 <- VariableFeaturePlot(pbmc.train) 
plot2 <- LabelPoints(plot = plot1, points = top10)
plot2 + guides(col="none")

## ----out.width="100%"---------------------------------------------------------
all.genes <- rownames(pbmc.train)
pbmc.train <- ScaleData(pbmc.train,features = all.genes)
pbmc.train <- RunPCA(pbmc.train, features = VariableFeatures(object = pbmc.train))
VizDimLoadings(pbmc.train, dims = 1:2, reduction = "pca")+theme(axis.text = element_text(size=10))
DimPlot(pbmc.train, reduction = "pca")+guides(col="none")

## -----------------------------------------------------------------------------
pbmc.train <- FindNeighbors(pbmc.train, dims = 1:10)
pbmc.train <- FindClusters(pbmc.train, resolution=0.5)

## -----------------------------------------------------------------------------
pbmc.train <- RunUMAP(pbmc.train, dims = 1:10)
DimPlot(pbmc.train, reduction = "umap")

## -----------------------------------------------------------------------------
clusters.train <- Idents(pbmc.train)
table(clusters.train)
length(clusters.train)
NCOL(Xtestsubset)

## -----------------------------------------------------------------------------
sf <- colSums(Xtestsubset)
Xtestsubset_norm <- t(apply(Xtestsubset, 1, function(u) u/sf))
Xtestsubset_lognorm <- log(Xtestsubset_norm +1)

## -----------------------------------------------------------------------------
cluster2 <- clusters.train==2
pvals2 <- apply(Xtestsubset_lognorm, 1, function(u) wilcox.test(u~cluster2)$p.value)
head(sort(pvals2))

## -----------------------------------------------------------------------------
pbmc.test <- pbmc.train
pbmc.test <- SetAssayData(object = pbmc.test, slot = "counts", new.data = Xtestsubset)
pbmc.test <- NormalizeData(pbmc.test)
pbmc.test <- ScaleData(pbmc.test, features = all.genes)

## -----------------------------------------------------------------------------
cluster2.markers <- FindMarkers(pbmc.test, ident.1 = 2, min.pct = 0)
head(sort(pvals2), n=10)
head(cluster2.markers, n = 10)

## -----------------------------------------------------------------------------
p0 <- DimPlot(pbmc.train, reduction = "umap")
p1 <- FeaturePlot(pbmc.train, features = c("RPS12","S100A9"), reduction="umap")
p2 <- FeaturePlot(pbmc.test, features = c("RPS12", "S100A9"), reduction="umap")
p0+p1+p2+plot_layout(nrow=3)

## -----------------------------------------------------------------------------
DimHeatmap(pbmc.train, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(pbmc.test, dims = 1:15, cells = 500, balanced = TRUE)

