#' Poisson count splitting
#' 
#' Takes one matrix of counts and returns a list containing two matrices of counts: a training set and a test set.
#'
#' The training set and the test set are independent under the assumption that the original data follow a Poisson distribution. 
#' Be sure to see the countsplitting tutorial vignette for details of how to use this correctly with existing
#' single cell RNA-seq pipelines.
#'
#' @export
#' @importFrom stats rbinom
#'
#' @param X A cell-by-gene matrix of integer counts
#' @param epsilon The thinning parameter for count splitting. Must be between 0 and 1.
#' @param folds The number of independent folds of data to return. Currently, the only option is 2. 
countsplit <- function(X, epsilon=0.5, folds=2) {
  if (epsilon <= 0 | epsilon >= 1) {
    stop('Parameter epsilon must be in (0,1)')
  }
  
  if (folds != 2) {
    stop("Splitting with more than two folds is coming soon.")
  }

  ## sparse formulation which samples only non-zero entries
  ## and replace sampled value inplace
  if (inherits(X,"dgCMatrix")) {
    Xtrain <- X
    Xtrain@x <- as.numeric(rbinom(n=length(X@x), size=X@x, prob=epsilon))
    Xtest <- X-Xtrain
  } else {
    ## dense formulation
    Xtrain <-apply(X,2,function(u) rbinom(n=length(u), size=u, prob=epsilon))
    rownames(Xtrain) <- rownames(X)
    colnames(Xtrain) <- colnames(X)
    Xtest <- X-Xtrain
    rownames(Xtest) <- rownames(X)
    colnames(Xtest) <- colnames(X)
  }
  return(list(train=Xtrain, test=Xtest))
}


#' @importFrom stats rbeta
betaBinSample <- function(x, b,eps) {
  p <- rbeta(length(x),eps*b,(1-eps)*b)
  return(rbinom(length(x),x,p))
}

#' Negative binomial count splitting
#' 
#' Takes one matrix of counts and returns a list containing two matrices of counts: a training set and a test set.
#'
#' The training set and the test set are independent under the assumption that the original data follow a Negative binomial distribution. 
#' Be sure to see the countsplit tutorial vignette for details of how to use this correctly with existing
#' This function is still in testing mode. Tutorials for this function are coming soon.
#'
#' @importFrom stats rbinom
#' @importFrom Matrix drop0
#' 
#' @export
#'
#' @param X A cell-by-gene matrix of integer counts
#' @param epsilon The thinning parameter for count splitting. Must be between 0 and 1.
#' @param overdisps A vector of length p, where p is the number of columns in X. This vector stores the 
#' gene specific overdispersion parameters. Note that a value of overdisp=infinity corresponds to the Poisson distribution. 
#' @param folds The number of independent folds of data to return. Currently, the only option is 2. 
nb.countsplit <- function(X, epsilon=0.5, overdisps=NULL, folds=2) {
  if (epsilon <= 0 | epsilon >= 1) {
    stop('Parameter epsilon must be in (0,1)')
  }
  
  ### If no overdispersion parameters are provided, just call the Poisson count split function.
  if (is.null(overdisps)) {
    return(countsplit(X, epsilon))
  }
  
  if (folds != 2) {
    stop("Splitting with more than two folds is coming soon.")
  }
  
  if (inherits(X,"dgCMatrix")){
    ##### Overdisp_ij = b_j if X_ij is nonzero, and is 0 otherwise.
    ## sparse formulation
    Xtrain <- X
    mapped_overdisps <- overdisps[Matrix::which(Xtrain != 0, arr.ind=T)[,"row"]] # for only the non-zero entries in X, maps the associated (gene-specific) overdispersion param
    Xtrain@x <- as.numeric(betaBinSample(Xtrain@x, mapped_overdisps, epsilon))
    Xtrain <- drop0(Xtrain)
    Xtest <- drop0(X - Xtrain)
    
  } else {
    ## dense formulation
    Xtrain <- t(apply(X, 1, function(u) betaBinSample(u, overdisps, epsilon)))
    Xtest <- X-Xtrain
    rownames(Xtrain) <- rownames(X)
    colnames(Xtrain) <- colnames(X)
    rownames(Xtest) <- rownames(X)
    colnames(Xtest) <- colnames(X)
  }
  return(list(train=Xtrain, test=Xtest))
}

