#' Takes one matrix of counts and returns two matrices of counts: a training set and a test set. 
#' 
#' Be sure to see the countsplitting tutorial vignette for details of how to use this correctly with existing 
#' single cell RNA-seq pipelines. 
#' 
#' @param X A cell-by-gene matrix of integer counts
#' @param epsilon The thinning parameter for count splitting. Must be between 0 and 1.
countsplit <- function(X, epsilon=0.5) {
  if (epsilon <= 0 | epsilon >= 1) {
    stop('Parameter epsilon must be in (0,1)')
  }
  
  Xtrain <-apply(X,2,function(u) rbinom(n=length(u), size=u, p=epsilon))
  rownames(Xtrain) <- rownames(X)
  colnames(Xtrain) <- colnames(X)
  Xtest <- X-Xtrain
  rownames(Xtest) <- rownames(X)
  colnames(Xtest) <- colnames(X)
  
  return(list(train=Xtrain, test=Xtest))
}


#' Anna note to self- should I make this called "testDE" or something and include t-tests and Wilcoxon tests as options??
#' IDK how much I should plan to do for the user. 
#' 
#' @importFrom speedglm speedglm
#'
#' @param counts A cell-by-gene matrix of integer counts.
#' @param latentvar The latent variable that you wish to test for differential expression across. Must have length equal to number of cells
#' @param family The family for the GLM. Typically will be "poisson", "quasipoisson", or "nb". Other options may cause errors.
#' @param test Specify the type of test you would like. For now the only choice is "wald". 
#' @param offsets Specify any offsets that you would like in your model. Typically a vector of logged size factors. Must have length equal to the number of cells. 
fit_glms <- function(counts, latentvar, family="poisson", test="wald", offsets=rep(1, NROW(counts))) {
  if (family=="nb") {
    pvals_pseudotime <- t(apply(Xtest, 2, function(u) mynbglm(u,clusters)))
   } 
  if (family=="quasipoisson") {
    res <- t(apply(counts, 2, function(u) as.numeric(summary(speedglm::speedglm(u~latentvar, family=quasipoisson("log")))$coefficients[2,])))
    }
  if (family=="poisson") {
    res <- t(apply(counts, 2, function(u) as.numeric(summary(speedglm::speedglm(u~latentvar, family=poisson("log")))$coefficients[2,])))
  }
  
  res <- data.frame(res[,c(1,2,4)])
  names(res) <- c("estimate", "se", "pval")
  return(res)
}

mynbglm <- function(u, latent, offset) {
  try1 <- try(suppressWarnings(MASS::glm.nb(u~latent, offset=offset)))
  if (class(try1)[1]=="try-error") {
    return(c(NA,NA,NA,NA))
  } else {
    return(as.numeric(summary(try1)$coefficients[2,]))
  }
}
