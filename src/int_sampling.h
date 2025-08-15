#ifndef INT_SAMPLING_H
#define INT_SAMPLING_H

#include <RcppEigen.h>

Rcpp::IntegerVector int_sampling(const int & K,
                                 const int & num_samples,
                                 const Rcpp::NumericVector & probs);

#endif  // INT_SAMPLING_H
