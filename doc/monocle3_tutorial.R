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

## -----------------------------------------------------------------------------
cds_train <- preprocess_cds(cds_train, num_dim = 50)
cds_train <- align_cds(cds_train, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

## -----------------------------------------------------------------------------
cds_train <- reduce_dimension(cds_train)
plot_cells(cds_train, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")

## -----------------------------------------------------------------------------
cds_train <- cluster_cells(cds_train)
cds_train <- learn_graph(cds_train)
plot_cells(cds_train,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

## -----------------------------------------------------------------------------
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds_train <- order_cells(cds_train, root_pr_nodes=get_earliest_principal_node(cds_train))
plot_cells(cds_train,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

## -----------------------------------------------------------------------------
cds_test <- cds_train
counts(cds_test) <- Xtest

## -----------------------------------------------------------------------------
ciliated_cds_pr_test_res <- graph_test(cds_test, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.0005))

## ---- results='hide'----------------------------------------------------------
plot_cells(cds_test, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

## -----------------------------------------------------------------------------
AFD_genes <- c("gcy-8", "dac-1", "oig-8")
AFD_lineage_cds <- cds_test[rowData(cds_test)$gene_short_name %in% AFD_genes,
                       colData(cds_test)$cell.type %in% c("AFD")]
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="embryo.time.bin",
                         min_expr=0.5)

