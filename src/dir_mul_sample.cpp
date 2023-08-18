#include <Rcpp.h>

// This function is based on code from the Rcpp Gallery https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/.
// It draws one realization from a multinomial distribution with size parameter "size" and probability vector "probs".
// Returns a single multinomial realization, whose length is equal to the number of elements in "probs".
// [[Rcpp::export]]
Rcpp::IntegerVector rmultinom_1(int &size, Rcpp::NumericVector &probs) {
  int N = probs.length();
  Rcpp::IntegerVector outcome(N);
  rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}

// Returns a single Dirichlet-Multinomial realization, who length is equal to "folds".
// The size parameter for the Dirichlet-Multinomial distribution is given by x.
// The other parameter vector has length folds, and each element is equal to b/folds.
// Assumes that you are allocating information equally between folds.
// Slower version gets called if that is not the case.
// [[Rcpp::export]]
Rcpp::IntegerVector dir_mul_sample_cpp(int &x, int folds, double b) {

  // Main goal is to intialize a vector p with length folds, which will serve as the input to multinomial sampling.

  // The easiest option is to have p=rep(1/folds, folds), which corresponds to regular multinomial sampling
  // AKA poisson count spliing. Initialize to this, and then update if b is not infinity.                                                          ./
  Rcpp::NumericVector ps(folds, 1.0/folds);

  // If b is 0, we want to randomly select one index of p to be 1, and all other indices should be 0
  // This is ugly C++ code, will consider updating.
  if (b==0) {
    // A strange case that shouldn't come up much.
    Rcpp::IntegerVector indices = Rcpp::seq_len(folds);
    int idx = Rcpp::as<int>(Rcpp::sample(indices, 1, false));

    for (int i = 0; i < ps.size(); i++) {
      ps[i] = 0;
    }
    ps[idx - 1] = 1;
  } else if (!(std::isinf(b))) {
    // This is the most important case, where we actually need to let p be a Dirichlet random variable
    // We generate the Dirichlet as a scaled gamma vector.
    Rcpp::NumericVector gammas = Rcpp::rgamma(folds, 1.0/folds * b, 1.0);
    double sum_gammas = std::accumulate(gammas.begin(), gammas.end(), 0.0);
    ps = gammas / sum_gammas;

    // Deal with this case, which can happen due to numerical error.
    // Revert to b=0 case if this happens.
    // I think this is the correct case to revert to, because it happens when b is very close to 0 but not exactly 0
    if (sum_gammas == 0) {
      Rcpp::IntegerVector indices = Rcpp::seq_len(folds);
      int idx = Rcpp::as<int>(Rcpp::sample(indices, 1, false));

      for (int i = 0; i < ps.size(); i++) {
        ps[i] = 0;
      }
      ps[idx - 1] = 1;
    }
  }

  /// why is this a matrix and not a vector? Check later
  Rcpp::IntegerVector result = rmultinom_1(x, ps);
  return result;
}

// So dumb it should jsut be done inline ofc.
Rcpp::IntegerVector mul_slower(int &x, Rcpp::NumericVector epsilon) {

}

// Literally just call something else if its Poisson, that will be faster.
Rcpp::IntegerVector dir_mul_slower(int &x, Rcpp::NumericVector epsilon, double b) {

  if (std::isinf(b)) {
    Rcpp::IntegerVector result = rmultinom_1(x, epsilon);
    return result;
  }

  // Main goal is to intialize a vector p with length folds, which will serve as the input to multinomial sampling.
  // The easiest option is to have p=rep(1/folds, folds), which corresponds to regular multinomial sampling
  // AKA poisson count spliing. Initialize to this, and then update if b is not infinity.                                                          ./
  int folds = epsilon.length();
  Rcpp::NumericVector ps(folds);

  // If b is 0, we want to randomly select one index of p to be 1, and all other indices should be 0
  // This is ugly C++ code, will consider updating.
  if (b==0) {
    // A strange case that shouldn't come up much.
    Rcpp::IntegerVector indices = Rcpp::seq_len(folds);
    int idx = Rcpp::as<int>(Rcpp::sample(indices, 1, false, epsilon));

    for (int i = 0; i < ps.size(); i++) {
      ps[i] = 0;
    }
    ps[idx - 1] = 1;
  } else {
    // This is the most important case, where we actually need to let p be a Dirichlet random variable
    // We generate the Dirichlet as a scaled gamma vector.
    Rcpp::NumericVector gammas(folds);
    for (int j=0; j<folds; j++) {
      gammas[j] = Rcpp::as<double>(Rcpp::rgamma(1, epsilon[j]*b, 1.0));
    }
    double sum_gammas = std::accumulate(gammas.begin(), gammas.end(), 0.0);
    ps = gammas / sum_gammas;

    // Deal with this case, which can happen due to numerical error.
    // Revert to b=0 case if this happens.
    // I think this is the correct case to revert to, because it happens when b is very close to 0 but not exactly 0
    if (sum_gammas == 0) {
      Rcpp::IntegerVector indices = Rcpp::seq_len(folds);
      int idx = Rcpp::as<int>(Rcpp::sample(indices, 1, false, epsilon));

      for (int i = 0; i < ps.size(); i++) {
        ps[i] = 0;
      }
      ps[idx - 1] = 1;
    }
  }

  /// why is this a matrix and not a vector? Check later
  Rcpp::IntegerVector result = rmultinom_1(x, ps);
  return result;
}


// Returns a single draw from a beta binomial distribution.
// This is the univariate version of the dirichlet multinomial, but implementing it as a separate case was
// useful to allow for a user-inputted epsilon.
// Returns a single draw from BetaBinomial(x, epsilon*b, (1-epsilon)*b)
// [[Rcpp::export]]
int beta_bin_sample_cpp(int &x, double eps, double b) {

  if (std::isinf(b)) {
    int x1 = Rcpp::as<int>(Rcpp::rbinom(1,x, eps));
    return(x1);
  }

  if (b==0) {
    Rcpp::IntegerVector choices = {x,0};
    Rcpp::NumericVector probs = {eps, 1-eps};
    int x1 = Rcpp::as<int>(Rcpp::sample(choices, 1, false, probs));
    return(x1);
  }

  Rcpp::NumericVector p = Rcpp::rbeta(1, eps * b, (1.0-eps)*b);
  double realp = Rcpp::as<double>(p);
  int x1 = Rcpp::as<int>(Rcpp::rbinom(1,x, realp));
  return(x1);
}


// Handles the efficient calling of dir_mul_sample across all elements in x
// [[Rcpp::export]]
Rcpp::IntegerMatrix mapply_dir_mul_sample_cpp(Rcpp::IntegerVector x, int folds, Rcpp::NumericVector overdisps) {
    int n = x.size();
    Rcpp::IntegerMatrix result(folds, n);
    Rcpp::IntegerVector sample(folds);

    for (int i = 0; i < n; i++) {
      sample = dir_mul_sample_cpp(x[i], folds, overdisps[i]);
      for (int j = 0; j < folds; j++) {
        result(j, i) = sample[j];
      }
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix mapply_dir_mul_slower(Rcpp::IntegerVector x, Rcpp::NumericVector epsilon, Rcpp::NumericVector overdisps) {
  int n = x.size();
  int folds = epsilon.length();
  Rcpp::IntegerMatrix result(folds, n);
  Rcpp::IntegerVector sample(folds);

  for (int i = 0; i < n; i++) {
    sample = dir_mul_slower(x[i], epsilon, overdisps[i]);
    for (int j = 0; j < folds; j++) {
      result(j, i) = sample[j];
    }
  }
  return result;
}

// Handles the efficient calling of beta_bino_sample across all elements in x
// [[Rcpp::export]]
Rcpp::IntegerMatrix mapply_betabin_sample_cpp(Rcpp::IntegerVector x, double eps1, Rcpp::NumericVector overdisps) {
  int n = x.size();
  Rcpp::IntegerMatrix result(2, n);
  Rcpp::IntegerVector sample(2);

  for (int i = 0; i < n; i++) {

    int x1 = beta_bin_sample_cpp(x[i], eps1, overdisps[i]);
    result(0, i) = x1;
    result(1,i) = x[i] - x1;

  }
  return result;
}
