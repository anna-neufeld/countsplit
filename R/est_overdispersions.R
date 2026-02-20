#' Estimate overdispersion parameters
#'
#' Takes one matrix of counts and estimates the overdispersion parameter of each row.
#'
#' 
#' @param matrix A genes (rows) x cells (columns) count matrix.
#' @param min_cells Minimum number of cells expressing a gene for it to be modeled by `sctransform::vst()`.
#' @param n_genes Number of genes to use in `vst()`. The default, `NULL`, uses all eligible genes.
#' @param ... Additional arguments passed to `sctransform::vst()`.
#'
#' @return A numeric vector of overdispersion values per row of matrix; genes not modeled return `Inf`.
#' 
#' @importFrom Matrix colSums rowSums
#' @importFrom stats setNames
#' @export
#' 
#' @examples
#' if (requireNamespace("sctransform", quietly = TRUE)) {
#'   library(Matrix)
#' n1 <- 500
#' n2 <- 500
#' p  <- 200
#' s1 <- 3     # True overdispersion parameter of first set of rows
#' s2 <- 9     # True overdispersion parameter of second set of rows
#' mean1 <- 9
#' mean2 <- 9
#' 
#'  Make both matrices with the  columns (same cells)
#' col_names <- paste0("cell_", seq_len(p))
#' X1 <- matrix(rnbinom(n1 * p, mu = mean1, size = s1), nrow = n1, ncol = p, dimnames = list(paste0("gene_", seq_len(n1)), col_names))
#' X2 <- matrix(rnbinom(n2 * p, mu = mean2, size = s2), nrow = n2, ncol = p, dimnames = list(paste0("gene_", (n1 + 1):(n1 + n2)), col_names))
#' 
#'  Stack by rows such that the resultant matrix is 1000 x 200
#' X <- rbind(X1, X2)
#' 
#' Estimate the overdispersions of the rows of the combined expression matrix
#' estover <- est_overdispersions(X)
#' 
#' Create a vector containing the true overdispersion parameters of first and second set of rows
#' v1 <- rep(s1, times = n1)
#' v2 <- rep(s2, times = n2)
#' trueover <- c(v1, v2)
#' 
#' Calculate the correlation between the estimated and true overdispersions
#' cor(estover, trueover)
#' countsplit(X, estover)


est_overdispersions <- function(matrix,
                                min_cells = 5,
                                n_genes = NULL, ...) {
  
  # Row names for output
  row_nms <- rownames(matrix)
  if (is.null(row_nms)) row_nms <- paste0("G", seq_len(nrow(matrix)))
  
  # Inf for everything as default (meaning "couldn't be calculated")
  out <- setNames(rep(Inf, nrow(matrix)), row_nms)
  
  # Drop all-zero columns (cells with zero total counts), since they cause issues w/ sctransform
  keep_cols <- Matrix::colSums(matrix) > 0
  if (!any(keep_cols)) return(out)
  matC <- matrix[, keep_cols, drop = FALSE] ### consider dropping this
  
  # Identify nonzero rows (genes with any counts)
  keep_rows <- Matrix::rowSums(matC) > 0
  if (!any(keep_rows)) return(out)
  mat2 <- matC[keep_rows, , drop = FALSE]
  
  # Run sctransform: genes are rows, cells are columns
  vst_res <- tryCatch(
    sctransform::vst(
      mat2,
      n_genes   = n_genes,
      min_cells = min_cells,
      ...
    ),
    error = function(e) NULL
  )
  
  # If vst failed, error message is returned
  if (is.null(vst_res) || is.null(vst_res$model_pars) || nrow(vst_res$model_pars) == 0) {
    message("CAUTION!: est_overdispersions: sctransform::vst() failed; returning Inf for all rows.")
    return(out)
  }
  
  # Extract overdispersion per gene
  theta <- as.numeric(vst_res$model_pars[, 1])
  theta_names <- rownames(vst_res$model_pars)
  
  # Any NA/NaN/Inf becomes Inf
  bad_idx <- !is.finite(theta)
  n_bad <- sum(bad_idx)
  
  if (n_bad > 0) {
    n_na  <- sum(is.na(theta))
    n_inf <- sum(is.infinite(theta))
    n_nan <- sum(is.nan(theta))  # NaN is also NA, but we report it explicitly
    message(
      "est_overdispersions: detected ", n_bad,
      " non-finite theta values (NA=", n_na,
      ", NaN=", n_nan,
      ", Inf=", n_inf,
      "); replacing with Inf."
    )
    theta[bad_idx] <- Inf
  }
  
  
  # Fill output for genes that got model_pars (others stay Inf)
  idx <- intersect(theta_names, names(out))
  out[idx] <- theta[match(idx, theta_names)]
  
  out
}

