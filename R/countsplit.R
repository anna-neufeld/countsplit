#' Takes one matrix of counts and returns a list containing two matrices of counts: a training set and a test set. 
#' 
#' Be sure to see the countsplitting tutorial vignette for details of how to use this correctly with existing 
#' single cell RNA-seq pipelines. 
#' 
#' @export
#' @importFrom stats rbinom
#' 
#' @param X A cell-by-gene matrix of integer counts
#' @param epsilon The thinning parameter for count splitting. Must be between 0 and 1.
countsplit <- function(X, epsilon=0.5) {
  if (epsilon <= 0 | epsilon >= 1) {
    stop('Parameter epsilon must be in (0,1)')
  }
  
  Xtrain <-apply(X,2,function(u) rbinom(n=length(u), size=u, prob=epsilon))
  rownames(Xtrain) <- rownames(X)
  colnames(Xtrain) <- colnames(X)
  Xtest <- X-Xtrain
  rownames(Xtest) <- rownames(X)
  colnames(Xtest) <- colnames(X)
  
  return(list(train=Xtrain, test=Xtest))
}

