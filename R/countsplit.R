#' Poisson count splitting
#' 
#' Takes one matrix of counts and returns a list containing two matrices of counts: a training set and a test set.
#'
#' The training set and the test set are independent under the assumption that the original data follow a Poisson distribution. 
#' Be sure to see the countsplitting tutorial vignette for details of how to use this correctly with existing
#' single cell RNA-seq pipelines.
#' 
#' To split your Poisson data matrix into multiple independent folds of data, please use the function "multisplit". 
#' 
#'
#' @export
#' @importFrom stats rbinom
#'
#' @param X A cell-by-gene matrix of integer counts
#' @param epsilon The thinning parameter for count splitting. Must be between 0 and 1.
countsplit <- function(X, folds=2, epsilon=rep(1/folds,folds),overdisps=rep(Inf, NCOL(X))) {
  return(multisplit(X, folds, epsilon, overdisps))
  # if (epsilon <= 0 | epsilon >= 1) {
  #   stop('Parameter epsilon must be in (0,1)')
  # }
  # 
  # ## sparse formulation which samples only non-zero entries
  # ## and replace sampled value inplace
  # if (inherits(X,"dgCMatrix")) {
  #   Xtrain <- X
  #   Xtrain@x <- as.numeric(rbinom(n=length(X@x), size=X@x, prob=epsilon))
  #   Xtest <- X-Xtrain
  # } else {
  #   ## dense formulation
  #   Xtrain <-apply(X,2,function(u) rbinom(n=length(u), size=u, prob=epsilon))
  #   rownames(Xtrain) <- rownames(X)
  #   colnames(Xtrain) <- colnames(X)
  #   Xtest <- X-Xtrain
  #   rownames(Xtest) <- rownames(X)
  #   colnames(Xtest) <- colnames(X)
  # }
  # return(list(train=Xtrain, test=Xtest))
}


#' Beta binomial sampling.
#'
#' The actual splitting funciton used to split a single negative binomial random variable into two independent pieces.
#' 
#" @importFrom stats rbeta
betaBinSample <- function(x, b,eps) {
  p <- rbeta(length(x),eps*b,(1-eps)*b)
  return(rbinom(length(x),x,p))
}

#' @importFrom stats rgamma
dirMulSample <- function(x, folds,b, epsilon=1/folds) {
  gammas <- rgamma(folds,epsilon*b,1)
  if (sum(gammas)==0) {
    gammas[sample(1:length(gammas),1)] <- 1
  }
  if (is.infinite(sum(gammas))) {
    gammas <- rep(1, length(gammas))
  }
  ps <- gammas/sum(gammas)
  x_partition <- rmultinom(1, x, ps)
  return(x_partition)
}

#' Negative binomial count splitting
#' 
#' Takes one matrix of counts and returns a list containing two matrices of counts: a training set and a test set.
#'
#' The training set and the test set are independent under the assumption that the original data follow a Negative binomial distribution,
#' with true overdispersion parameter equal to the "overdisps" argument. Be sure to see the countsplit tutorial vignette for details of how to use this correctly with existing scRNA-seq packages. 
#'
#' To split your negative binomial data matrix into multiple independent folds of data, please use the function "multisplit". 
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
nb.countsplit <- function(X, epsilon=0.5, overdisps=NULL) {
  return(multisplit(X, folds, epsilon, overdisps))
  
  #if (epsilon <= 0 | epsilon >= 1) {
  #  stop('Parameter epsilon must be in (0,1)')
  #}
  
  #if (is.null(overdisps)) {
  #  message("Since no overdispersion parameters were provided, we are performing Poisson count splitting.")
  #  return(countsplit(X, epsilon, folds))
  #}
  
  #if (folds != 2) {
  #  stop("Splitting with more than two folds is coming soon.")
  #}
  
  #if (inherits(X,"dgCMatrix")){
    ##### Overdisp_ij = b_j if X_ij is nonzero, and is 0 otherwise.
    ## sparse formulation
   # Xtrain <- X
  #  mapped_overdisps <- overdisps[Matrix::which(Xtrain != 0, arr.ind=T)[,"row"]] # for only the non-zero entries in X, maps the associated (gene-specific) overdispersion param
  #  Xtrain@x <- as.numeric(betaBinSample(Xtrain@x, mapped_overdisps, epsilon))
   # Xtrain <- drop0(Xtrain)
  #  Xtest <- drop0(X - Xtrain)
    
  #} else {
  #  ## dense formulation
  #  Xtrain <- t(apply(X, 1, function(u) betaBinSample(u, overdisps, epsilon)))
  #  Xtest <- X-Xtrain
  #  rownames(Xtrain) <- rownames(X)
  #  colnames(Xtrain) <- colnames(X)
  #  rownames(Xtest) <- rownames(X)
  #  colnames(Xtest) <- colnames(X)
  #}
  #return(list(train=Xtrain, test=Xtest))
}




#' Multifold count splitting
#' 
#' Takes one matrix of counts and returns a list containing several matrices of counts, which can be used to do cross validation.
#' If the parameter `overdisps` is set to NULL, this performs multifold Poisson count splitting.
#' If the parameter `overdisps` is not null, then this perfoms multifold negative binomial count splitting.
#' Please see our tutorial website for details on how to use this function properly. 
#' 
#' Returns a three dimensional matrix, which has dimensions folds by cells by genes. Or maybe a list of sparse matrices. 
#'
#' @export
#' @importFrom stats rbinom
#'
#' @param X A cell-by-gene matrix of integer counts. For the moment, must be a dense matrix. The sparse implementation is coming soon. 
#' @param folds The number of folds of data to return.
#' @param epsilon A vector whose length is equal to folds. The components all must be between 0 and 1 and must sum to 1. The most common (and the detault) option is to have every entry of epsilon set to 1/folds. 
#' @param overdisps If NULL of Inf, then this function assumes that you would like to do Poisson count splitting. Otherwise, overdisps should be a vector with length equal to the number of columns of X. The entries should speccify the gene-specific overdispersion parameters to be used for splitting.
multisplit <- function(X, folds=2, epsilon=rep(1/folds, folds), overdisps = NULL) {
  if ( (sum(epsilon) != 1) | (sum(epsilon < 0) != 0)) {
    stop('The elements of the epsilon vector must be non-zero and must sum to 1.')
  }
  
  if (is.null(overdisps)) {
    overdisps <- rep(Inf, ncol(X))
  }
  
  ### Do everything as sparse matrices, for speed. 
  if (!(inherits(X,"dgCMatrix"))) {
    X <-  as(X, "sparseMatrix") 
  }
   
  mapped_overdisps <- overdisps[Matrix::which(X  != 0, arr.ind=T)[,"col"]] # for only the non-zero entries in X, maps the associated (gene-specific) overdispersion param
  results <- mapply(function(x,y) as.numeric(dirMulSample(x, folds, y, epsilon)), X@x, mapped_overdisps)
    
  partition <- vector(mode = "list", length = folds)
  for (f in 1:folds) {
    Xfold <- X
    Xfold@x <- as.numeric(results[f,])
    partition[[f]] <- Xfold
  }
  
  return(partition)
}
