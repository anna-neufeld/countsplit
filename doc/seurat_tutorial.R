## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("Seurat")

## -----------------------------------------------------------------------------
remotes::install_github("anna-neufeld/countsplit")

## -----------------------------------------------------------------------------
library(Seurat)
library(countsplit)
library(ggplot2)
library(patchwork)
library(mclust)
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
pbmc <- CreateSeuratObject(counts = pbmc.counts, min.cells = 3, min.features = 200)

## -----------------------------------------------------------------------------
# Apply to training object
pbmc.train[["percent.mt"]] <- PercentageFeatureSet(pbmc.train, pattern = "^MT-")
pbmc.train <- subset(pbmc.train, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Apply to full object
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

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
dim(pbmc.train)
dim(pbmc)

## -----------------------------------------------------------------------------
all.equal(GetAssayData(pbmc.train, "counts"), GetAssayData(pbmc.train, "data"))

## -----------------------------------------------------------------------------
GetAssayData(pbmc.train, "scale.data")

## -----------------------------------------------------------------------------
pbmc.train <- NormalizeData(pbmc.train)
all.equal(GetAssayData(pbmc.train, "counts"), GetAssayData(pbmc.train, "data"))

# Apply to pbmc, for the sake of comparison. 
pbmc <- NormalizeData(pbmc)

## -----------------------------------------------------------------------------
dim(pbmc.train)
pbmc.train <- FindVariableFeatures(pbmc.train, selection.method = "vst", nfeatures = 2000)
dim(pbmc.train)

# Apply to pbmc, for the sake of comparison. 
pbmc  <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

## -----------------------------------------------------------------------------
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc) + ggtitle("pbmc")
plot2 <- LabelPoints(plot = plot1, points = top10)
top10.train <- head(VariableFeatures(pbmc.train), 10)
plot1.train <- VariableFeaturePlot(pbmc.train) + ggtitle("pbmc.train")
plot2.train <- LabelPoints(plot = plot1.train, points = top10.train)


plot2 + plot2.train & guides(col="none")



## ----out.width="100%"---------------------------------------------------------
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
p1 <- VizDimLoadings(pbmc, dims = 1, reduction = "pca")+theme(axis.text = element_text(size=10))+ggtitle("pbmc")
p2 <- VizDimLoadings(pbmc, dims = 2, reduction = "pca")+theme(axis.text = element_text(size=10))

all.genes.train <- rownames(pbmc.train)
pbmc.train <- ScaleData(pbmc.train,features = all.genes.train)
pbmc.train <- RunPCA(pbmc.train, features = VariableFeatures(object = pbmc.train))
p1.train <- VizDimLoadings(pbmc.train, dims = 1, reduction = "pca")+theme(axis.text = element_text(size=10))+ggtitle("pbmc.train")
p2.train <- VizDimLoadings(pbmc.train, dims = 2, reduction = "pca")+theme(axis.text = element_text(size=10))

p1+p1.train+p2+p2.train+plot_layout(nrow=2, ncol=2)

## -----------------------------------------------------------------------------
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution=0.5)
pbmc.train <- FindNeighbors(pbmc.train, dims = 1:10)
pbmc.train <- FindClusters(pbmc.train, resolution=0.5)

## -----------------------------------------------------------------------------
pbmc.train <- RunUMAP(pbmc.train, dims = 1:10)
DimPlot(pbmc.train, reduction = "umap")

## -----------------------------------------------------------------------------
clusters.train <- Idents(pbmc.train)
clusters.full.subset <- Idents(pbmc)[colnames(pbmc.train)]
length(clusters.train)
length(clusters.full.subset)
adjustedRandIndex(clusters.train, clusters.full.subset)

## -----------------------------------------------------------------------------
table(clusters.train, clusters.full.subset)

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
DimHeatmap(pbmc.train, dims = 1:15, cells = 500, balanced = TRUE)

## -----------------------------------------------------------------------------
DimHeatmap(pbmc.test, dims = 1:15, cells = 500, balanced = TRUE)

