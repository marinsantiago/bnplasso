// -----------------------------------------------------------------------------
// Update mu
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "update_mu.h"


/**
 * Update the intercept term, mu.                                              \
 *                                                                             \
 * @param X Eigen::MatrixXd: X matrix.                                         \
 * @param y Eigen::VectorXd: y vector.                                         \
 * @param beta Eigen::VectorXd: Current value of the beta parameters.          \
 * @param sigma2 double: Current value of sigma2.                              \
 *                                                                             \
 * @return double: Updated value of mu.                                        \
 */
// [[Rcpp::export]]
double update_mu(const Eigen::MatrixXd & X,
                 const Eigen::VectorXd & y,
                 const Eigen::VectorXd & beta,
                 const double sigma2) {
  
  const double n = y.size();
  const double eps = 1e-08;
  const Eigen::VectorXd res = y - X * beta;
  const double post_mean = res.array().mean();
  const double post_sd = std::sqrt(std::max(sigma2 / n, eps));
  return R::rnorm(post_mean, post_sd);
}
