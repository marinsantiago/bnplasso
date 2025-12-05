// -----------------------------------------------------------------------------
// Utils
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "utils.h"


/**
 * Function to convert an Rcpp::NumericVector into an Eigen::VectorXd          \ 
 *                                                                             \
 * @param x Rcpp::NumericVector: A vector to convert into Eigen::VectorXd      \
 *                                                                             \
 * @return Eigen::VectorXd: Converted vector.                                  \
 */
// [[Rcpp::export]]
Eigen::VectorXd NumVec_to_EigenVec(Rcpp::NumericVector & x) {
  Eigen::Map<Eigen::VectorXd> out(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(x));
  return out;
}


/**
 * Function to convert an Rcpp::NumericVector into an Eigen::VectorXf (float)! \ 
 *                                                                             \
 * @param x Rcpp::NumericVector: A vector to convert into Eigen::VectorXf      \
 *                                                                             \
 * @return Eigen::VectorXf: Converted vector (as float)!                       \
 */
// [[Rcpp::export]]
Eigen::VectorXf NumVec_to_EigenVec_float(const Rcpp::NumericVector & x) {
  Eigen::VectorXf out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    out[i] = static_cast<float>(x[i]);
  }
  return out;
}


/**
 * Function to convert an Eigen::VectorXd into an Rcpp::NumericVector          \ 
 *                                                                             \
 * @param x Eigen::VectorXd: A vector to convert into Rcpp::NumericVector      \
 *                                                                             \
 * @return Rcpp::NumericVector: Converted vector.                              \
 */
// [[Rcpp::export]]
Rcpp::NumericVector EigenVec_to_NumVec(const Eigen::VectorXd & x) {
  Rcpp::NumericVector out(x.data(), x.data() + x.size());
  return out;
}


/**
 *                                                                             \
 * Function to compute the cross-product t(X) %*% X                            \ 
 *                                                                             \
 */
// [[Rcpp::export]]
Eigen::MatrixXd get_tXX(const Eigen::MatrixXd & X) {
  Eigen::MatrixXd out = X.transpose() * X;
  return out;
}
