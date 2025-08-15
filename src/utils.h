#ifndef UTILS_H
#define UTILS_H

#include <RcppEigen.h>

Eigen::VectorXd NumVec_to_EigenVec(Rcpp::NumericVector & x);

Rcpp::NumericVector EigenVec_to_NumVec(const Eigen::VectorXd & x);

Eigen::MatrixXd get_tXX(const Eigen::MatrixXd & X);

#endif  // UTILS_H
