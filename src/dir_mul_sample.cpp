#include <Rcpp.h>

// This function is based on code from the Rcpp Gallery https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/.
// It draws one realization from a multinomial distribution with size parameter "size" and probability vector "probs".
// Returns a single multinomial realization, whose length is equal to the number of elements in "probs".
// @param size the multinomial size parameter. Should be a positive integer.
// @param probs the multinomial probability parameter. Should be a vector with non-negative entries that sum to 1.
// [[Rcpp::export]]
Rcpp::IntegerVector rmultinom_1(int &size, Rcpp::NumericVector &probs) {
  int N = probs.length();
  Rcpp::IntegerVector outcome(N);
  rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}

// Returns a single Dirichlet-Multinomial realization, who length is equal to "folds".
// Assumes that we want all folds to have equal information in them- slower version must be called if this is not the case.
// The size parameter is "size", and the other parameters are all given by "b/folds".
// Corresponds to drawing a vector called "probs" from a Dirichlet distribution with all parameters equal to b/folds,
// and then drawing from a multinomial(size, probs) distribution.
// [[Rcpp::export]]
Rcpp::IntegerVector dir_mul_sample_cpp(int &size, int folds, double b) {

  // Sampling from a dirichlet multinomial distribution is the same as sampling from a multinomial distribution
  // with a randomly generated "probs" vector. The main work of this function involves generating the
  // appropriate "probs" vector, which has length "folds".

  // The easiest option is to have p=rep(1/folds, folds), which corresponds to regular multinomial ampling
  // with equally distributed folds. We will initialize to this, and then update if we don't
  // AKA poisson count splitting. Initialize to this, and then update if b is not infinity.                                                          ./
  Rcpp::NumericVector probs(folds, 1.0/folds);

  // When b is not infinity, need to sample "probs" from a Dirichlet distribution,
  // which can be accomplished by sampling from a gamma distribution and then normalizing.
  if (!std::isinf(b)) {

      Rcpp::NumericVector gammas = Rcpp::rgamma(folds, 1.0/folds * b, 1.0);
      double sum_gammas = std::accumulate(gammas.begin(), gammas.end(), 0.0);
      probs = gammas / sum_gammas;

      // The above code will fail if sum_gammas = 0.
      // This can happen in the case where b=0 or is very close to 0.
      // In this case, we sample from the limiting case of the Dirichlet distribution,
      // which sets one element of probs to 1 and sets all of the rest to 0.
      if (sum_gammas == 0) {
        // Is this efficient in C++??
        Rcpp::IntegerVector indices = Rcpp::seq_len(folds);
        int idx = Rcpp::as<int>(Rcpp::sample(indices, 1, false));
        for (int i = 0; i < probs.size(); i++) {
          probs[i] = 0;
        }
        probs[idx - 1] = 1;
    }
  }

  // Finally ready to sample from multinomial distribution and return the result.
  Rcpp::IntegerVector result = rmultinom_1(size, probs);
  return result;
  }


// Returns a single Dirichlet-Multinomial realization, who length is equal to "folds".
// Does not assume that all folds have equal information, so is slower than "dir_mul_sample_cpp".
// The size parameter is "size", and the other parameters are given by the vector "epsilon*b".
// Corresponds to drawing a vector called "probs" from a Dirichlet distribution with parameter "epsilon*b",
// and then drawing from a multinomial(size, probs) distribution.
Rcpp::IntegerVector dir_mul_slower(int &x, Rcpp::NumericVector epsilon, double b) {

  // In the case where b=infinity, we can be efficient.
  if (std::isinf(b)) {
    Rcpp::IntegerVector result = rmultinom_1(x, epsilon);
    return result;
  }

  // Otherwise, main goal is to initialize a vector probs with length folds, which will serve as the input to multinomial sampling.
  int folds = epsilon.length();
  Rcpp::NumericVector probs(folds);

  Rcpp::NumericVector gammas(folds);
  // The slow thing is that we must call rgamma once per fold.
  for (int j=0; j<folds; j++) {
    gammas[j] = Rcpp::as<double>(Rcpp::rgamma(1, epsilon[j]*b, 1.0));
  }
  double sum_gammas = std::accumulate(gammas.begin(), gammas.end(), 0.0);
  probs = gammas / sum_gammas;

  // Special case needed if sum_gammas = 0.
  // Do the limiting case.
  if (sum_gammas == 0) {
    Rcpp::IntegerVector indices = Rcpp::seq_len(folds);
    int idx = Rcpp::as<int>(Rcpp::sample(indices, 1, false, epsilon));
    for (int i = 0; i < probs.size(); i++) {
        probs[i] = 0;
    }
    probs[idx - 1] = 1;
    }

  // Finally ready to sample from multinomial.
  Rcpp::IntegerVector result = rmultinom_1(x, probs);
  return result;
}

// Returns a single draw from a beta binomial distribution.
// This is the univariate version of the Dirichlet Multinomial,
// but implementing it separately allows for additional efficiency.
// Returns a single draw from BetaBinomial(x, epsilon*b, (1-epsilon)*b)
// [[Rcpp::export]]
int beta_bin_sample_cpp(int &x, double eps, double b) {

  // Special case.
  if (std::isinf(b)) {
    int x1 = Rcpp::as<int>(Rcpp::rbinom(1,x, eps));
    return(x1);
  }

  // Special case
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
// Entry j in x gets turned into a vector with length folds, where these entries
// were drawn from a Dirichlet-multinomial distribution with parameters x and overdisps[j]/folds.
// Assumes that we are allocating information equally between folds.
// If this is not the case, the slower version must be called.
// [[Rcpp::export]]
Rcpp::IntegerMatrix mapply_dir_mul_sample_cpp(Rcpp::IntegerVector x, int folds, Rcpp::NumericVector overdisps) {
    int n = x.size();
    Rcpp::IntegerMatrix result(folds, n);
    Rcpp::IntegerVector sample(folds);

    // Is this seriously the most efficient way to do this?
    for (int i = 0; i < n; i++) {
      sample = dir_mul_sample_cpp(x[i], folds, overdisps[i]);
      for (int j = 0; j < folds; j++) {
        result(j, i) = sample[j];
      }
    }
    return result;
}

// Handles the efficient calling of dir_mul_sample across all elements in x.
// The length of the vector epsilon tells us how many folds we are generating.
// Entry j in x gets turned into a vector with length folds, where these entries
// were drawn from a Dirichlet-multinomial distribution with parameters x and epsilon*overdisps[j].
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

// Handles the efficient calling of beta_bin_sample_cpp across all elements in x.
// Entry j in x gets turned into (x1,x2), where x1 is drawn from a betabinomial distribution with parameters
// x, eps1*overdisp[j], and (1-eps)*overdisps[j]. And x2=x-x1.
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
