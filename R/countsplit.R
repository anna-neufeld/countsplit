#' Count splitting
#'
#' Takes one matrix of counts and splits it into a specified number of folds. Each fold is a matrix of counts with the same dimension
#' as the original matrix. Summing element-wise across the folds yields the original data matrix.
#'
#' When the argument `overdisps` is set to NULL, this function performs the Poisson count splitting methodology outlined in
#' Neufeld et al. (2022). With this setting, the folds of data are independent only if the original data were drawn from a Poisson distribution.
#'
#' If the data are thought to be overdispersed relative to the Poisson, then we may instead model them as coming from a negative binomial distribution,
#' If we assume that \eqn{X_{ij} \sim NB(\mu_{ij}, b_j)}, where this parameterization means that \eqn{ E[X_{ij}] = \mu_{ij}} and \eqn{ Var[X_{ij}] = \mu_{ij} + \mu_{ij}^2/b_j}, then
#' we should pass in `overdisps` = \eqn{c(b_1, \ldots, b_j)}. If this is the correct assumption, then the resulting folds of data will be independent.
#' This is the negative binomial count splitting method of Neufeld et al. (2023).
#'
#' Please see our tutorials and vignettes for more details.
#'
#' @export
#' @importFrom methods as
#' @import Matrix
#' @import Rcpp
#' @useDynLib countsplit
#' @param X A cell-by-gene matrix of integer counts. Note that this differs from many scRNA-seq packages, where gene-by-cell is the convention. When Poisson count splitting is used (`overdisps=NULL`), the matrix can be either cell-by-gene or gene-by-cell.
#' @param folds An integer specifying how many folds you would like to split your data into.
#' @param epsilon A vector, which has length `folds`, that stores non-zero elements that sum to one. Determines the proportion of information from X that is allocated to each fold.
#' When `folds` is not equal to 2, the recommended (and default) setting is to allocate equal amounts of information to each fold, such that each element is `1/folds`.
#' When `folds=2`, the default is still `(1/2, 1/2)`, but other values may be beneficial.
#' @param overdisps If NULL, then Poisson count splitting will be performed. Otherwise, this parameter should be a vector of non-negative numbers whose length is equal to the number of columns of `X`.
#' These numbers are the overdispersion parameters for each column in `X`. If these are unknown, they can be estimated with a function such as
#' `vst` in the package `sctransform`.
#' @return A list of length `folds`. Each element in the list stores a sparse matrix with the same dimensions as the data `X`. Each list element is a fold of data.
#'
#' @examples
#' library(countsplit)
#' library(Matrix)
#' library(Rcpp)
#' # A Poisson count splitting example.
#' n=400
#' p=2
#' X <- matrix(rpois(n*p, 7), nrow=n, ncol=p)
#' split <- countsplit(X, folds=2)
#' Xtrain <- split[[1]]
#' Xtest <- split[[2]]
#' cor(Xtrain[,1], Xtest[,1])
#' cor(Xtrain[,2], Xtest[,2])
#'
#' # A negative binomial count splitting example.
#' X <- matrix(rnbinom(n*p, mu=7, size=7), nrow=n, ncol=p)
#' split <- countsplit(X, folds=2, overdisps=c(7,7))
#' Xtrain <- split[[1]]
#' Xtest <- split[[2]]
#' cor(Xtrain[,1], Xtest[,1])
#' cor(Xtrain[,2], Xtest[,2])
#' @references reference
countsplit <- function(X, folds=2, epsilon=rep(1/folds, folds), overdisps = NULL) {
  if (is.null(overdisps)) {
    overdisps <- rep(Inf, ncol(X))
    message("As no overdispersion parameters were provided, Poisson count splitting will be performed.")
  }

  if (length(overdisps) != NCOL(X)) {
    stop("You should provide one overdispersion parameter for every column of X. Make sure that your matrix X is cell-by-gene rather than gene-by-cell.")
  }
  
  
  if (length(epsilon) != folds | sum(epsilon) != 1 | sum(epsilon <= 0) > 0) {
    stop("The parameter epsilon should be a vector with length folds with positive entries that sum to 1.")
  }

  ### Do everything as sparse matrices, for speed.
  if (!(inherits(X,"dgCMatrix"))) {
    X <-  as(X, "sparseMatrix")
  }

  # For only the non-zero entries in X, maps the associated (gene-specific) overdispersion parameter so that the
  # entries correspond to X@x.
  mapped_overdisps <- overdisps[Matrix::which(X  != 0, arr.ind=T)[,"col"]]


  # Things are less efficient in the case where information is not equally allocated between folds.
  if (!identical(epsilon, rep(1/folds, folds))) {

    if (folds==2) {
      # Can still be efficient if only two folds.
      results <- mapply_betabin_sample_cpp(x = X@x, eps1 = epsilon[1], overdisps = mapped_overdisps) # use Rcpp function to increase the performance
    } else {
      # This is the only remaining slow case, and would be an unusual choice of parameters.
      results <- mapply_dir_mul_slower(x = X@x, epsilon = epsilon, overdisps = mapped_overdisps)
    }
  } else {
    # faster case where everything is equally allocated across folds.
    results <- mapply_dir_mul_sample_cpp(x = X@x, folds = folds, overdisps = mapped_overdisps) # use Rcpp function to increase the performance
  }

  # Reformat "results" into a list where each list element stores a fold of data.
  # Each fold of data is represented as a sparse matrix.
  partition <- vector(mode = "list", length = folds)
  for (f in 1:folds) {
    Xfold <- X
    Xfold@x <- as.numeric(results[f,])
    partition[[f]] <- Xfold
  }
  return(partition)
}


