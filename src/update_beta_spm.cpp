// -----------------------------------------------------------------------------
// Update beta within a sparse means problem under:
//   (a) An independent prior on sigma2.
//   (b) A conjugate prior on sigma2.
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "update_beta_spm.h"


/**
 * Function to generate random draws from a multivariate normal distribution,  \
 * using the algorithm proposed by                                             \
 * Rue (2001) <https://doi.org/10.1111/1467-9868.00288>, which is based on a   \
 * Cholesky factorization of the precision matrix.                             \
 * This function assumes an independent variance prior as in                   \
 * Moran et al. (2019) <doi:10.1214/19-BA1149>.                                \
 *                                                                             \
 * @param y Eigen::VectorXd: Vector of responses of size \eqn{n}.              \
 * @param tau2 Eigen::VectorXd: Current value of tau2.                         \
 * @param sigma2 double: Sampling variance.                                    \
 *                                                                             \
 * @return Eigen::VectorXd: Vector of regression coefficients drawn from its   \
 * full conditional distribution.                                              \
 */
// [[Rcpp::export]]
Eigen::VectorXd update_beta_spm_ind_sigma(const Eigen::VectorXd & y,
                                          const Eigen::VectorXd & tau2,
                                          const double sigma2) {
  
  // Pre-compute constants and prepare returns ---------------------------------
  const Eigen::Index p = tau2.size();
  const double sigma2_inv = 1.0 / sigma2;
  const double eps = 1e-15;
  Eigen::VectorXd out(p);
  
  // Iterate over entries ------------------------------------------------------
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    const double tau2_inv_j = 1.0 / tau2[jj];
    const double phi_inv_j = std::max(1.0 / (sigma2_inv + tau2_inv_j), eps);
    // Compute the mean
    const double m = phi_inv_j * y[jj] / sigma2; // In spm: n == p.
    // Compute the standard deviation
    const double sd = std::sqrt(phi_inv_j);
    // Generate a random draw
    out[jj] = R::rnorm(m, sd);
  }
  
  return out;
}


/**
 * Function to generate random draws from a multivariate normal distribution,  \
 * using the algorithm proposed by                                             \
 * Rue (2001) <https://doi.org/10.1111/1467-9868.00288>, which is based on a   \
 * Cholesky factorization of the precision matrix.                             \
 * This function assumes a conjugate prior on the variance. See                \
 * Moran et al. (2019) <doi:10.1214/19-BA1149> for details.                    \
 *                                                                             \
 * @param y Eigen::VectorXd: Vector of responses of size \eqn{n}.              \
 * @param tau2 Eigen::VectorXd: Current value of tau2.                         \
 * @param sigma2 double: Sampling variance.                                    \
 *                                                                             \
 * @return Eigen::VectorXd: Vector of regression coefficients drawn from its   \
 * full conditional distribution.                                              \
 */
// [[Rcpp::export]]
Eigen::VectorXd update_beta_spm_conj_sigma(const Eigen::VectorXd & y,
                                           const Eigen::VectorXd & tau2,
                                           const double sigma2) {
  
  // Pre-compute constants and prepare returns ---------------------------------
  const Eigen::Index p = tau2.size();
  const double eps = 1e-15;
  Eigen::VectorXd out(p);
  
  // Iterate over entries ------------------------------------------------------
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    const double tau2_inv_j = 1.0 / tau2[jj];
    const double phi_inv_j = 1.0 / (1.0 + tau2_inv_j);
    // Compute the mean
    const double m = phi_inv_j * y[jj]; // In spm: n == p.
    // Compute the standard deviation
    const double sd = std::sqrt(std::max(sigma2 * phi_inv_j, eps));
    // Generate a random draw
    out[jj] = R::rnorm(m, sd);
  }
  return out;
}
