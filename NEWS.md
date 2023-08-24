# countsplit 4.0.0

Thanks to some Github contributions from Mischko Heming, this version of the package has been substantially sped up compared to version 3.0.0. This speedup added a dependency on Rcpp, as the main function is now implemented in C++. 

# countsplit 3.0.0

This is the initial CRAN submission for countsplit. Compared to previous versions we have:

1. Moved the data and tutorials to a new package to save on build time and package size.
2. Fully implemented Poisson and negative binomial count splitting, both in their 2-fold and their multi-fold versions.
3. Streamlined so that all types of count splitting are included in the same `countsplit()` function. 
