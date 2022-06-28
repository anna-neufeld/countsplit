#' Cardiomyocyte data. 
#'
#' This is a subset of a dataset that was collected by Elborany et al., 2022. This data was also analyzed in Neufeld et al., 2022. 
#' This subset contains 10,000 cells and more than 20,000 genes. The format is a SingleCellExperiment object.
#'
#' @docType data
#'}
#' @keywords datasets
#' @format A RDS object. 
#' @examples
#' @references Elorbany, Reem, et al. "Single-cell sequencing reveals lineage-specific dynamic genetic regulation of gene expression during human cardiomyocyte differentiation." PLoS genetics 18.1 (2022): e1009666.
#' data(cm)
"cm"

#' PBMC data
#' 
#' This dataset is downloaded at the start of the Seurat guided clustering tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).
#' It is a dataset of Peripheral Blood Mononuclear Cells (PBMC) which has been made freely available by 10X Genomics. 
#' The dataset contains raw counts for 2,700 single cells that were sequenced on the Illumina NextSeq 500. 
#'
#' @docType data
#'}
#' @keywords datasets
#' @format A matrix. 
#' @examples
#' data(pbmc.counts)
"pbmc.counts"
