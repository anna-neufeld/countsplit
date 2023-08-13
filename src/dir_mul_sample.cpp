#include <Rcpp.h>
using namespace Rcpp;

IntegerVector rmultinom_1(int &size, NumericVector &probs, int &N) {
    IntegerVector outcome(N);
    rmultinom(size, probs.begin(), N, outcome.begin());
    return outcome;
}

IntegerMatrix rmultinom_rcpp(int n, int size, NumericVector &probs) {
    int N = probs.length();
    IntegerMatrix sim(N, n);
    for (int i = 0; i < n; i++) {
        sim(_,i) = rmultinom_1(size, probs, N);
    }
    return sim;
}

// [[Rcpp::export]]
IntegerMatrix dir_mul_sample_cpp(int &x, int folds, double b) {
  double epsilon = 1.0 / folds;
  NumericVector gammas = rgamma(folds, epsilon * b, 1.0);
  double sum_gammas = std::accumulate(gammas.begin(), gammas.end(), 0.0);

  if (sum_gammas == 0) {
        IntegerVector indices = seq_len(folds);
        int idx = Rcpp::as<int>(Rcpp::sample(indices, 1, false));
        gammas[idx - 1] = 1;
    }

  if (std::isinf(sum_gammas)) {
    gammas = NumericVector(gammas.length(), 1);
  }
  NumericVector ps = gammas / sum_gammas;
  IntegerMatrix result = rmultinom_rcpp(1, x, ps);

  return result;
}

// [[Rcpp::export]]
IntegerMatrix mapply_dir_mul_sample_cpp(IntegerVector x, int folds, NumericVector overdisps) {
    int n = x.size();
    IntegerMatrix result(folds, n);
    IntegerMatrix sample(folds, 1);

    for (int i = 0; i < n; i++) {
      sample = dir_mul_sample_cpp(x[i], folds, overdisps[i]);
      for (int j = 0; j < folds; j++) {
        result(j, i) = sample(j, 0);
      }
    }

    return result;
}
