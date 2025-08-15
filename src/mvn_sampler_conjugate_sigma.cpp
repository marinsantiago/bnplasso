// -----------------------------------------------------------------------------
// Multivariate normal samplers, assuming a conjugate prior on the variance
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include "utils.h"
#include "mvn_sampler_conjugate_sigma.h"


/**
 * Function to generate random draws from a multivariate normal distribution,  \
 * using the algorithm proposed by                                             \
 * Rue (2001) <https://doi.org/10.1111/1467-9868.00288>, which is based on a   \
 * Cholesky factorization of the precision matrix.                             \
 * This function assumes a conjugate prior on the variance. See                \
 * Moran et al. (2019) <doi:10.1214/19-BA1149> for details.                    \
 *                                                                             \
 * @param X Eigen::MatrixXd: Matrix of predictors of size \eqn{n}-by-\eqn{d},  \
 *  where each  of the \eqn{n} rows is an observation vector.                  \
 * @param tX Eigen::MatrixXd: Transpose of the matrix X.                       \
 * @param tXX Eigen::MatrixXd: Cross-product of the form: t(X) %*% X           \
 * @param y Eigen::VectorXd: Vector of responses of size \eqn{n}.              \
 * @param tau2 Eigen::VectorXd: Vector of hyper-parameters of the form:        \
 *  (tau_1^2, ..., tau_p^2)'.                                                  \
 * @param sigma2 double: Sampling variance.                                    \
 * @param mu double: Intercept term.                                           \
 *                                                                             \
 * @return Eigen::VectorXd: Vector of regression coefficients drawn from its   \
 *  full conditional distribution.
 */
// [[Rcpp::export]]
Eigen::VectorXd sample_beta_conj_sigma(const Eigen::MatrixXd & X,
                                       const Eigen::MatrixXd & tX,
                                       const Eigen::MatrixXd & tXX,
                                       const Eigen::VectorXd & y,
                                       const Eigen::VectorXd & tau2,
                                       const double sigma2,
                                       const double mu) {
  
  // Pre-compute constants and prepare the returns
  const int p = tau2.size();
  Eigen::VectorXd out;
  
  // Compute (tau_1^2, ..., tau_p^2)^(-1)
  Eigen::VectorXd tau2_inv_pos = tau2.array().inverse();
  tau2_inv_pos = tau2_inv_pos.array().abs() + 1e-08;
  // Compute Phi
  Eigen::MatrixXd tau2_inv_mat = tau2_inv_pos.asDiagonal();
  Eigen::MatrixXd Phi = tXX + tau2_inv_mat;
  // Compute X^t * (y - mu)
  Eigen::VectorXd y_mu = y.array() - mu;
  Eigen::VectorXd tXy = tX * y_mu;
  // Cholesky factorization
  Eigen::LLT<Eigen::MatrixXd> llt(Phi);
  // Solve the inner and outer systems using the Cholesky factorization
  Eigen::VectorXd mu_tilde = llt.solve(tXy);
  // Sample z ~ MVN(0, I_p)
  Rcpp::NumericVector z_Rcpp = Rcpp::rnorm(p, 0.0, 1.0);
  Eigen::VectorXd z = NumVec_to_EigenVec(z_Rcpp);
  // Sample beta
  Eigen::MatrixXd Lt = llt.matrixU();
  Eigen::VectorXd sys_solve = Lt.triangularView<Eigen::Upper>().solve(z);
  double sigma = std::sqrt(std::abs(sigma2) + 1e-8);
  sys_solve *= sigma;
  out = mu_tilde + sys_solve;

  return out;
}
