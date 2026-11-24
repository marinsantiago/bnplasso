// -----------------------------------------------------------------------------
// Update lambda under a Bayesian Lasso and a Bayesian adaptive Lasso
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "update_lambda_lasso.h"


/**
 * Update the lambda parameters under a Bayesian Lasso.                        \
 *                                                                             \
 * @param tau2 Eigen::VectorXd: Current value of the tau2 parameters.          \
 * @param a double: Shape parameter of the prior centering measure.            \
 * @param b double: Rate parameter of the prior centering measure.             \
 *                                                                             \
 * @return Eigen::VectorXd: Updated values of lambda2.                         \
 */
// [[Rcpp::export]]
Eigen::VectorXd update_lambda2_Bayesian_lasso(
    const Eigen::VectorXd & tau2,
    const double a, 
    const double b
) {
  const int p = tau2.size();
  const double eps = 1e-08;
  const double sum_tau2_half = 0.5 * tau2.sum();
  const double post_shape = std::max(a + p, eps);
  const double post_rate = std::max(b + sum_tau2_half, eps);
  const double draw = std::max(R::rgamma(post_shape, 1 / post_rate), eps);
  Eigen::VectorXd out(p);
  for (int jj = 0; jj < p; ++jj) { out[jj] = draw; }
  return out;
} 


/**
 * Update the lambda parameters under a Bayesian adaptive Lasso.               \
 *                                                                             \
 * @param tau2 Eigen::VectorXd: Current value of the tau2 parameters.          \
 * @param a double: Shape parameter of the prior centering measure.            \
 * @param b double: Rate parameter of the prior centering measure.             \
 *                                                                             \
 * @return Eigen::VectorXd: Updated values of lambda2.                         \
 */
// [[Rcpp::export]]
Eigen::VectorXd update_lambda2_Bayesian_adapt_lasso(
    const Eigen::VectorXd & tau2,
    const double a,
    const double b
) {
  const Eigen::Index p = tau2.size();
  const double eps = 1e-08;
  const double a_one = std::max(a + 1, eps);
  Eigen::VectorXd out(p);
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    const double post_shape = a_one;
    const double post_rate = std::max(b + 0.5 * tau2[jj], eps);
    out[jj] = std::max(R::rgamma(post_shape, 1 / post_rate), eps);
  }
  return out;
}
