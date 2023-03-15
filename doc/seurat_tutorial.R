## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("Seurat")

## ---- eval=FALSE--------------------------------------------------------------
#  remotes::install_github("anna-neufeld/countsplit")

## -----------------------------------------------------------------------------
library(Seurat)
library(countsplit)
library(ggplot2)
library(patchwork)
library(mclust)

## -----------------------------------------------------------------------------
data(pbmc.counts, package="countsplit")

## -----------------------------------------------------------------------------
rownames(pbmc.counts) <- sapply(rownames(pbmc.counts), function(u) stringr::str_replace_all(u, "_","-"))

## ---- collapse=TRUE-----------------------------------------------------------
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
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## ---- collapse=TRUE-----------------------------------------------------------
dim(Xtrain)
dim(Xtest)
dim(pbmc.train)

## ---- collapse=TRUE-----------------------------------------------------------
rows <- rownames(pbmc.train)
cols <- colnames(pbmc.train)
Xtestsubset <- Xtest[rows,cols]
dim(Xtestsubset)

## ---- collapse=TRUE-----------------------------------------------------------
dim(pbmc.train)
dim(pbmc)

## ---- collapse=TRUE-----------------------------------------------------------
all.equal(GetAssayData(pbmc.train, "counts"), GetAssayData(pbmc.train, "data"))
GetAssayData(pbmc.train, "scale.data")

## -----------------------------------------------------------------------------
pbmc.train <- NormalizeData(pbmc.train)
all.equal(GetAssayData(pbmc.train, "counts"), GetAssayData(pbmc.train, "data"))

## ---- collapse=TRUE, results='hide'-------------------------------------------
dim(pbmc.train)
pbmc.train <- FindVariableFeatures(pbmc.train, selection.method = "vst", nfeatures = 2000)
dim(pbmc.train)

## ---- results='hide'----------------------------------------------------------
all.genes.train <- rownames(pbmc.train)
pbmc.train <- ScaleData(pbmc.train,features = all.genes.train)
pbmc.train <- RunPCA(pbmc.train, features = VariableFeatures(object = pbmc.train))

## ---- results='hide'----------------------------------------------------------
pbmc <- NormalizeData(pbmc)
pbmc  <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

## -----------------------------------------------------------------------------
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc) + ggtitle("pbmc")
plot2 <- LabelPoints(plot = plot1, points = top10)
top10.train <- head(VariableFeatures(pbmc.train), 10)
plot1.train <- VariableFeaturePlot(pbmc.train) + ggtitle("pbmc.train")
plot2.train <- LabelPoints(plot = plot1.train, points = top10.train)
plot2 + plot2.train & guides(col="none")

## ---- collapse=TRUE-----------------------------------------------------------
sort(top10)
sort(top10.train)

## ----fig.height=9, fig.width=5------------------------------------------------
p1 <- VizDimLoadings(pbmc, dims = 1, reduction = "pca")+theme(axis.text = element_text(size=7))+ggtitle("pbmc")
p2 <- VizDimLoadings(pbmc, dims = 2, reduction = "pca")+theme(axis.text = element_text(size=7))
p1.train <- VizDimLoadings(pbmc.train, dims = 1, reduction = "pca")+theme(axis.text = element_text(size=7))+ggtitle("pbmc.train")
p2.train <- VizDimLoadings(pbmc.train, dims = 2, reduction = "pca")+theme(axis.text = element_text(size=7))

p1+p1.train+p2+p2.train+plot_layout(nrow=2, ncol=2)

## ---- results='hide'----------------------------------------------------------
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution=0.5)
pbmc.train <- FindNeighbors(pbmc.train, dims = 1:10)
pbmc.train <- FindClusters(pbmc.train, resolution=0.5)

## ---- results='hide'----------------------------------------------------------
pbmc.train <- RunUMAP(pbmc.train, dims = 1:10)
DimPlot(pbmc.train, reduction = "umap")

## -----------------------------------------------------------------------------
clusters.train <- Idents(pbmc.train)
clusters.full <- Idents(pbmc)
length(clusters.train)
length(clusters.full)

## -----------------------------------------------------------------------------
clusters.full.subset <- clusters.full[colnames(pbmc.train)] 
adjustedRandIndex(clusters.train, clusters.full.subset)

## -----------------------------------------------------------------------------
table(clusters.train, clusters.full.subset)

## -----------------------------------------------------------------------------
clusters.train <- Idents(pbmc.train)
length(clusters.train)
NCOL(Xtestsubset)

## -----------------------------------------------------------------------------
## Log normalize the test set
sf.test <- colSums(Xtestsubset)
Xtestsubset_norm <- t(apply(Xtestsubset, 1, function(u) u/sf.test))
Xtestsubset_lognorm <- log(Xtestsubset_norm +1)

## Log normalize the full dataset
Xsubset <- pbmc.counts[rownames(pbmc),colnames(pbmc)]
sf.full <- colSums(Xsubset)
Xsubset_norm <- t(apply(Xsubset, 1, function(u) u/sf.full))
Xsubset_lognorm <- log(Xsubset_norm +1)

## ---- collapse=TRUE-----------------------------------------------------------
### Do the count splitting analysis
cluster2.train <- clusters.train==2
pvals2.countsplit <- apply(Xtestsubset_lognorm, 1, function(u) wilcox.test(u~cluster2.train)$p.value)

### Do the naive method analysis 
clusters.full <- Idents(pbmc)
cluster2.full <- clusters.full==2
pvals2.naive <- apply(Xsubset_lognorm, 1, function(u) wilcox.test(u~cluster2.full)$p.value)

head(sort(pvals2.countsplit))
head(sort(pvals2.naive))

## ---- results='hide'----------------------------------------------------------
pbmc.train[['test']] <- CreateAssayObject(counts=Xtestsubset)
pbmc.train <- NormalizeData(pbmc.train, assay="test")
pbmc.train <- ScaleData(pbmc.train, assay="test")

## ---- collapse=TRUE-----------------------------------------------------------
cluster2.markers <- FindMarkers(pbmc.train, ident.1=2, min.pct=0, assay="test")

head(sort(pvals2.countsplit), n=10)
head(cluster2.markers, n = 10)

## ---- collapse=TRUE-----------------------------------------------------------
head(sort(pvals2.naive), n=10)
head(FindMarkers(pbmc, ident.1=2), n=10)

## ----fig.height=9, fig.width=6------------------------------------------------
DimHeatmap(pbmc.train, dims = 1:15, cells = 500, balanced = TRUE, nfeatures=10)

## ---- fig.height=9, fig.width=6-----------------------------------------------
DimHeatmap(pbmc.train, dims = 1:15, cells = 500, balanced = TRUE, nfeatures=10, assay="test")

