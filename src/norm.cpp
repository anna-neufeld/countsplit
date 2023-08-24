#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]



//'  modified Log normalization
//'
//' This is log normalization function based on Seurats normalization method and slightly altered for this package.
//' Reference: https://github.com/satijalab/seurat/blob/763259d05991d40721dee99c9919ec6d4491d15e/src/data_manipulation.cpp#L113
//'
//' @param data A sparse gene-by-cell matrix of integer counts
//' @param scale_factor a scale factor (usually 10000)
//' @param eps_train A double that determines the proportion of information that is allocated to each fold. Default is 1 - 1/folds.
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> LogNorm(Eigen::SparseMatrix<double> data, int scale_factor, double eps_train){
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  for (int k=0; k < data.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      it.valueRef() = std::log1p((double(it.value())/ eps_train) / colSums[k] * scale_factor);
    }
  }
  return data;
}
