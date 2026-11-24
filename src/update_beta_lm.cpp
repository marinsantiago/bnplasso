// -----------------------------------------------------------------------------
// Update beta within a linear model under:
//   (a) An independent prior on sigma2.
//   (b) A conjugate prior on sigma2.
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "update_beta_lm.h"


/**
 * Function to generate random draws from a multivariate normal distribution,  \
 * using the algorithm proposed by                                             \
 * Rue (2001) <https://doi.org/10.1111/1467-9868.00288>, which is based on a   \
 * Cholesky factorization of the precision matrix.                             \
 * This function assumes an independent variance prior as in                   \
 * Moran et al. (2019) <doi:10.1214/19-BA1149>.                                \
 *                                                                             \
 * @param Xt Eigen::MatrixXd: Transpose of the matrix X.                       \
 * @param XtX Eigen::MatrixXd: Cross-product of the form: t(X) %*% X           \
 * @param y Eigen::VectorXd: Vector of responses of size \eqn{n}.              \
 * @param tau2 Eigen::VectorXd: Current value of tau2.                         \
 * @param sigma2 double: Sampling variance.                                    \
 * @param mu double: Intercept term.                                           \
 *                                                                             \
 * @return Eigen::VectorXd: Vector of regression coefficients drawn from its   \
 * full conditional distribution.                                              \
 */
// [[Rcpp::export]]
Eigen::VectorXd update_beta_lm_ind_sigma(const Eigen::MatrixXd & Xt,
                                         const Eigen::MatrixXd & XtX,
                                         const Eigen::VectorXd & y,
                                         const Eigen::VectorXd & tau2,
                                         const double sigma2,
                                         const double mu) {
  
  // Pre-compute constants -----------------------------------------------------
  const Eigen::Index p = tau2.size();
  const double eps = 1e-08;
  
  // Compute phi matrix --------------------------------------------------------
  const Eigen::VectorXd tau2_inv = tau2.array().inverse();
  Eigen::MatrixXd phi = XtX / sigma2;
  phi += tau2_inv.asDiagonal();
  phi.diagonal().array() += eps;
  
  // Compute mean vector -------------------------------------------------------
  const Eigen::VectorXd y_mu = y.array() - mu;
  const Eigen::VectorXd mean_vec = Xt * y_mu / sigma2;
  
  // Cholesky factorization ----------------------------------------------------
  Eigen::LLT<Eigen::MatrixXd> llt(phi);
  // Solve the inner and outer systems using the Cholesky factorization
  const Eigen::VectorXd mu_tilde = llt.solve(mean_vec);
  
  // Sample z ~ MVN(0, I_p) ----------------------------------------------------
  Eigen::VectorXd z(p);
  for (Eigen::Index jj = 0; jj < p; ++jj) { z[jj] = R::rnorm(0.0, 1.0); }
  
  // Sample beta ---------------------------------------------------------------
  const Eigen::MatrixXd Lt = llt.matrixU();
  return mu_tilde + Lt.triangularView<Eigen::Upper>().solve(z);
}


/**
 * Function to generate random draws from a multivariate normal distribution,  \
 * using the algorithm proposed by                                             \
 * Rue (2001) <https://doi.org/10.1111/1467-9868.00288>, which is based on a   \
 * Cholesky factorization of the precision matrix.                             \
 * This function assumes a conjugate prior on the variance. See                \
 * Moran et al. (2019) <doi:10.1214/19-BA1149> for details.                    \
 *                                                                             \
 * @param tX Eigen::MatrixXd: Transpose of the matrix X.                       \
 * @param tXX Eigen::MatrixXd: Cross-product of the form: t(X) %*% X           \
 * @param y Eigen::VectorXd: Vector of responses of size \eqn{n}.              \
 * @param tau2 Eigen::VectorXd: Current value of tau2.                         \
 * @param sigma2 double: Sampling variance.                                    \
 * @param mu double: Intercept term.                                           \
 *                                                                             \
 * @return Eigen::VectorXd: Vector of regression coefficients drawn from its   \
 * full conditional distribution.                                              \
 */
// [[Rcpp::export]]
Eigen::VectorXd update_beta_lm_conj_sigma(const Eigen::MatrixXd & Xt,
                                          const Eigen::MatrixXd & XtX,
                                          const Eigen::VectorXd & y,
                                          const Eigen::VectorXd & tau2,
                                          const double sigma2,
                                          const double mu) {
  
  // Pre-compute constants -----------------------------------------------------
  const Eigen::Index p = tau2.size();
  const double eps = 1e-08;
  
  // Compute phi matrix --------------------------------------------------------
  const Eigen::VectorXd tau2_inv_pos = tau2.array().inverse();
  Eigen::MatrixXd phi = XtX;
  phi += tau2_inv_pos.asDiagonal();
  phi.diagonal().array() += eps;  

  // Compute mean vector -------------------------------------------------------
  const Eigen::VectorXd y_mu = y.array() - mu;
  const Eigen::VectorXd mean_vec = Xt * y_mu;
  
  // Cholesky factorization ----------------------------------------------------
  Eigen::LLT<Eigen::MatrixXd> llt(phi);
  // Solve the inner and outer systems using the Cholesky factorization
  const Eigen::VectorXd mu_tilde = llt.solve(mean_vec);
  
  // Sample z ~ MVN(0, I_p) ----------------------------------------------------
  Eigen::VectorXd z(p);
  for (Eigen::Index jj = 0; jj < p; ++jj) { z[jj] = R::rnorm(0.0, 1.0); }
  
  // Sample beta ---------------------------------------------------------------
  const double sigma = std::sqrt(std::max(sigma2, eps));
  const Eigen::MatrixXd Lt = llt.matrixU();
  Eigen::VectorXd sys_solve = Lt.triangularView<Eigen::Upper>().solve(z);
  sys_solve *= sigma;
  return mu_tilde + sys_solve;
}
