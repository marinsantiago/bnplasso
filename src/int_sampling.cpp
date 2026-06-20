// -----------------------------------------------------------------------------
// Random (integer) sampling with replacement
// -----------------------------------------------------------------------------

#include <Rcpp.h>

// [[Rcpp::depends(Rcpp)]]
using namespace Rcpp;

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

/**
 * Random sampling with replacement.                                           \
 *                                                                             \
 * @param K int: A positive number, the number of items to choose from.        \
 * @param num_samples int: A non-negative integer giving the number of items   \
 *  to choose.                                                                 \
 * @param probs NumericVector: a vector of probability weights for obtaining   \
 *  the elements of the vector being sampled.                                  \
 *                                                                             \
 * @return IntegerVector: An integer vector of length `num_samples` with       \
 *  elements from 1:K                                                          \
 */
// [[Rcpp::export]]
Rcpp::IntegerVector int_sampling(const int & K,
                                 const int & num_samples,
                                 const Rcpp::NumericVector & probs) {
  
  // Prepare the returns
  Rcpp::IntegerVector out(num_samples);
  
  // Pre-compute constants
  Rcpp::NumericVector cum_probs = Rcpp::cumsum(probs);
  Rcpp::NumericVector random_values = Rcpp::runif(num_samples);
  auto cum_probs_begin = cum_probs.begin();
  auto cum_probs_end = cum_probs.end();
  
  for (int i = 0; i < num_samples; ++i) {
    out[i] = std::lower_bound(cum_probs_begin,
                              cum_probs_end,
                              random_values[i]) - cum_probs_begin + 1;
  }
  
  return out;
}
