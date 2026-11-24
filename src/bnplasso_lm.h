#ifndef BNPLASSO_LM_H
#define BNPLASSO_LM_H

#include <RcppEigen.h>

Rcpp::List bnplasso_lm(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const bool intercept,
    Rcpp::Nullable<Rcpp::String> shrinkage_prior,
    Rcpp::Nullable<Rcpp::NumericVector> a,
    Rcpp::Nullable<Rcpp::NumericVector> b,
    Rcpp::Nullable<Rcpp::NumericVector> alpha,
    Rcpp::Nullable<Rcpp::String> variance_prior_type,
    const int max_iters,
    const int burn_in,
    const int thin,
    const bool polya_urn_sampler,
    const bool verbose
);

#endif  // BNPLASSO_LM_H
