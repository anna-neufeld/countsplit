#include <Rcpp.h>

// the first two function were based on the Rcpp Gallery https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/
// to mimic the R function rmultinom
Rcpp::IntegerVector rmultinom_1(int &size, Rcpp::NumericVector &probs, int &N) {
    Rcpp::IntegerVector outcome(N);
    rmultinom(size, probs.begin(), N, outcome.begin());
    return outcome;
}

Rcpp::IntegerMatrix rmultinom_rcpp(int n, int size, Rcpp::NumericVector &probs) {
  int N = probs.length();
  Rcpp::IntegerMatrix sim(N, n);
  for (int i = 0; i < n; i++) {
        sim(Rcpp::_,i) = rmultinom_1(size, probs, N);
    }
    return sim;
  }

// [[Rcpp::export]]
Rcpp::IntegerMatrix dir_mul_sample_cpp(int &x, int folds, double b) {
  double epsilon = 1.0 / folds;
  Rcpp::NumericVector gammas = Rcpp::rgamma(folds, epsilon * b, 1.0);
  double sum_gammas = std::accumulate(gammas.begin(), gammas.end(), 0.0);

  if (sum_gammas == 0) {
        Rcpp::IntegerVector indices = Rcpp::seq_len(folds);
        int idx = Rcpp::as<int>(Rcpp::sample(indices, 1, false));
        gammas[idx - 1] = 1;
    }

  if (std::isinf(sum_gammas)) {
    gammas = Rcpp::NumericVector(gammas.length(), 1);
  }
  Rcpp::NumericVector ps = gammas / sum_gammas;
  Rcpp::IntegerMatrix result = rmultinom_rcpp(1, x, ps);

  return result;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix mapply_dir_mul_sample_cpp(Rcpp::IntegerVector x, int folds, Rcpp::NumericVector overdisps) {
    int n = x.size();
    Rcpp::IntegerMatrix result(folds, n);
    Rcpp::IntegerMatrix sample(folds, 1);

    for (int i = 0; i < n; i++) {
      sample = dir_mul_sample_cpp(x[i], folds, overdisps[i]);
      for (int j = 0; j < folds; j++) {
        result(j, i) = sample(j, 0);
      }
    }

    return result;
}
