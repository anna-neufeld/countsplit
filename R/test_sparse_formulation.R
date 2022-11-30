library(countsplit)

## from https://anna-neufeld.github.io/countsplit/articles/seurat_tutorial.html
data(pbmc.counts, package="countsplit")
rownames(pbmc.counts) <- sapply(rownames(pbmc.counts), function(u) stringr::str_replace_all(u, "_","-"))

#########################################
## 1) run countsplit w/ sparse formulation
startTime <- Sys.time()
set.seed(1)
split <- countsplit(pbmc.counts, epsilon=0.5)
Xtrain <- split$train
Xtest <- split$test

endTime <- Sys.time()
print(endTime - startTime) #Time difference of 0.3524549 secs


#########################################
## 2) run countsplit w/ dense formulation
startTime2 <- Sys.time()
pbmc.counts.dense = as.matrix(pbmc.counts)
set.seed(1)
split2 <- countsplit(pbmc.counts.dense, epsilon=0.5)
Xtrain.dense <- split2$train
Xtest.dense <- split2$test

endTime2 <- Sys.time()
print(endTime2 - startTime2) # Time difference of 6.527817 secs
