#ifndef UPDATE_LAMBDA_BNP_H
#define UPDATE_LAMBDA_BNP_H

#include <RcppEigen.h>

Rcpp::List update_lambda_DP_Polya_urn(
    const Eigen::VectorXd & tau2,
    const Eigen::VectorXd & lambda2_k,
    const Rcpp::IntegerVector & clust_indicators,
    const double a, 
    const double b,
    const double alpha
);

Rcpp::List update_lambda_DP_blocked_Gibbs(
    const Eigen::VectorXd & tau2,
    const Eigen::VectorXd & lambda2_candidates,   
    const Rcpp::IntegerVector & clust_indicators, 
    const Eigen::VectorXd & omega,
    const double a,
    const double b,
    const double alpha
);

#endif  // UPDATE_LAMBDA_BNP_H
