// -----------------------------------------------------------------------------
// Random integer sampling
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "int_sampling.h"


/**
 * Random sampling with replacement.                                           \
 *                                                                             \
 * @param K int: A positive number, the number of items to choose from.        \
 * @param num_samples int: A non-negative integer giving the number of items   \
 *  to choose.                                                                 \
 * @param probs NumericVector: a vector of probabilities (or un-normalized     \ 
 *  weights) for obtaining the elements of the vector being sampled.           \
 *                                                                             \
 * @return IntegerVector: An integer vector of length `num_samples` with       \
 *  elements from 1:K                                                          \
 */
// [[Rcpp::export]]
Rcpp::IntegerVector int_sampling(const int K,
                                 const int num_samples,
                                 const Eigen::VectorXd & probs) {
  
  // Prepare the returns -------------------------------------------------------
  Rcpp::IntegerVector out(num_samples);
  
  // Build the cumulative distribution -----------------------------------------
  Rcpp::NumericVector cum_probs(K);
  double* const cum_pointer = cum_probs.begin(); 
  double total = 0.0;
  for (int k = 0; k < K; ++k) {
    const double p = probs[k];
    total += p;
    cum_pointer[k] = total;
  }
  
  // Integer sampling ----------------------------------------------------------
  int* const out_pointer = out.begin();
  double* const cum_end = cum_pointer + K;
  for (int i = 0; i < num_samples; ++i) {
    const double u = R::unif_rand() * total;
    const std::ptrdiff_t idx = std::lower_bound(
      cum_pointer, cum_end, u
    ) - cum_pointer; 
    out_pointer[i] = static_cast<int>(idx) + 1;
  }

  return out;
}
