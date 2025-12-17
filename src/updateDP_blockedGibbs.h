#ifndef UPDATEDP_BLOCKEDGIBBS_H
#define UPDATEDP_BLOCKEDGIBBS_H

#include <RcppEigen.h>

Rcpp::IntegerVector step2_blockedGibbs_dbl(const Eigen::VectorXd & tau2,
                                           const Eigen::VectorXd & lambda2_all,
                                           const Eigen::VectorXd & omega,
                                           const int K);

Rcpp::IntegerVector step2_blockedGibbs_flt(const Eigen::VectorXf & tau2,
                                           const Eigen::VectorXf & lambda2_all,
                                           const Eigen::VectorXf & omega,
                                           const int K);

#endif  // UPDATEDP_BLOCKEDGIBBS_H
