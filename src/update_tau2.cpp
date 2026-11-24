// -----------------------------------------------------------------------------
// Update tau2 within a linear model under:
//   (a) An independent prior on sigma2.
//   (b) A conjugate prior on sigma2.
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "rinvgauss.h"
#include "update_tau2.h"


/**
 * Update tau2 within a linear model under an independent prior on sigma2.     \
 *                                                                             \
 * @param beta Eigen::VectorXd: Current value of the beta parameters.          \
 * @param lambda2 Eigen::VectorXd: Current value of the lambda2 parameters.    \
 *                                                                             \
 * @return Eigen::VectorXd: Updated values of tau2.                            \
 */
// [[Rcpp::export]]
Eigen::VectorXd update_tau2_ind_sigma(const Eigen::VectorXd & beta,
                                      const Eigen::VectorXd & lambda2) {
  
  const Eigen::Index p = beta.size();
  const double eps = 1e-08;
  Eigen::VectorXd out(p);
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    const double l2 = lambda2[jj];
    const double b2 = beta[jj] * beta[jj];
    const double xi = std::sqrt(std::max(l2 / b2, eps));
    const double invGaussian = rinvgauss(xi, l2);
    out[jj] = std::max(1.0 / invGaussian, eps);
  }
  return out;
}


/**
 * Update tau2 within a linear model under a conjugate prior on sigma2.        \
 *                                                                             \
 * @param beta Eigen::VectorXd: Current value of the beta parameters.          \
 * @param lambda2 Eigen::VectorXd: Current value of the lambda2 parameters.    \
 * @param sigma2 double: Current value of sigma2s.                             \
 *                                                                             \
 * @return Eigen::VectorXd: Updated values of tau2.                            \
 */
// [[Rcpp::export]]
Eigen::VectorXd update_tau2_conj_sigma(const Eigen::VectorXd & beta,
                                       const Eigen::VectorXd & lambda2,
                                       const double sigma2) {
  
  const Eigen::Index p = beta.size();
  const double eps = 1e-08;
  Eigen::VectorXd out(p);
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    const double l2 = lambda2[jj];
    const double l2_s2 = l2 * sigma2;
    const double b2 = beta[jj] * beta[jj];
    const double xi = std::sqrt(std::max(l2_s2 / b2, eps));
    const double invGaussian = rinvgauss(xi, l2);
    out[jj] = std::max(1.0 / invGaussian, eps);
  }
  return out;
}
