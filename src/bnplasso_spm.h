#ifndef BNPLASSO_SPM_H
#define BNPLASSO_SPM_H

#include <RcppEigen.h>

Rcpp::List bnplasso_spm(
    const Eigen::VectorXd & y,
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

#endif  // BNPLASSO_SPM_H
