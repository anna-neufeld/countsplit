#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' @name LogNorm
//' @title modified Log normalization
//' @description This function performs log(x+1) normalization after divided by the total counts per cell and multiplied by a scale factor, as done in Seurat.
//' The only difference is that we here multiply by 1/eps_train.
//' Also this function has been rewritten in RcppArmadillo.
//' @param mat A sparse gene-by-cell matrix of integer counts
//' @param scale_factor a scale factor (usually 10000)
//' @param eps_train A double that determines the proportion of information that is allocated to each fold. Default is 1 - 1/folds.
//' @return A sparse matrix of log-normalized counts multiplied by 1/eps_train.
//' @export
// [[Rcpp::export]]
arma::sp_mat LogNorm(arma::sp_mat mat, int scale_factor, double eps_train) {
  arma::uword numCols = mat.n_cols;
  arma::vec colSums(numCols);

  for(arma::sp_mat::const_iterator it = mat.begin(); it != mat.end(); ++it) {
    colSums(it.col()) += *it;
  }

  for (arma::sp_mat::iterator it = mat.begin(); it != mat.end(); ++it) {
    *it = std::log1p((*it / eps_train) / colSums[it.col()] * scale_factor);
  }
  return mat;
}
