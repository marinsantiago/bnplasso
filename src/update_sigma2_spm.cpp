// -----------------------------------------------------------------------------
// Update sigma2 within a sparse means problem under:
//   (a) An independent prior on sigma2.
//   (b) A conjugate prior on sigma2.
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "update_sigma2_spm.h"


/**
 * Sample from the full conditional distribution of sigma2, assuming an        \
 * independent variance prior.                                                 \
 */
// [[Rcpp::export]]
double update_sigma2_spm_ind_prior(const Eigen::VectorXd & y,
                                   const Eigen::VectorXd & beta) {

  const int n = y.size();
  const double eps = 1e-08;
  const Eigen::VectorXd res = y - beta;
  const double ssr = res.transpose() * res;
  const double a_post = std::max(0.5 * n, eps);
  const double b_post = std::max(0.5 * ssr, eps);
  return 1.0 / std::max(R::rgamma(a_post, 1.0 / b_post), eps);
}


/**
 * Sample from the full conditional distribution of sigma2, assuming a         \
 * conjugate variance prior.                                                   \
 */
// [[Rcpp::export]]
double update_sigma2_spm_conj_prior(const Eigen::VectorXd & y,
                                    const Eigen::VectorXd & beta,
                                    const Eigen::VectorXd & tau2) {
  const double n = y.size();
  const double p = beta.size();
  const double eps = 1e-08;
  // Quadratic form of betas
  Eigen::VectorXd inv_tau2 = tau2.array().inverse();
  inv_tau2.array() += eps;
  const Eigen::VectorXd beta_inv_tau2 = beta.array() * inv_tau2.array();
  const double quad_form_beta = beta_inv_tau2.transpose() * beta;
  // Residuals and inverse-gamma parameters
  const Eigen::VectorXd res = y - beta;
  const double ssr = res.transpose() * res;
  const double a_post = std::max(0.5 * (n + p), eps);
  const double b_post = std::max(0.5 * (ssr + quad_form_beta), eps);
  return 1.0 / std::max(R::rgamma(a_post, 1.0 / b_post), eps);
}
