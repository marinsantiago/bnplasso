// -----------------------------------------------------------------------------
// Multivariate normal samplers, assuming a conjugate prior on the variance
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include "updateDP_blockedGibbs.h"
#include "int_sampling.h"
#include "utils.h"


/**
 * Function to perform "Step 2" in the blocked Gibbs update of the underlying  \ 
 * Dirichlet process mixture in the nonparametric Bayesian Lasso.              \
 *                                                                             \
 * @param tau2 VectorXd: Vector of size \eqn{p} with current values for tau2.  \
 * @param lambda2_all VectorXd: Vector of size \eqn{K} with all the current    \
 *  possible values for lambda2                                                \
 * @param omega VectorXd:VectorXd: Vector of size \eqn{K} with all the current \
 *  probabilities from each component.                                         \
 * @param K int: Truncation.                                                   \
 *                                                                             \
 * @return Rcpp::IntegerVector: Vector of updated cluster allocations          \
 */
// [[Rcpp::export]]
Rcpp::IntegerVector step2_blockedGibbs_dbl(const Eigen::VectorXd & tau2,
                                           const Eigen::VectorXd & lambda2_all,
                                           const Eigen::VectorXd & omega,
                                           const int K) {
  
  // Pre-compute constants and prepare the returns
  const int p = tau2.size();
  Rcpp::IntegerVector out(p);
  
  for (int j = 0; j < p; ++j) {
    Eigen::VectorXd prob_j = omega; // Initialize vector of probabilities
    for (int k = 0; k < K; ++k) {
      double l_2k = lambda2_all[k] * 0.5;
      prob_j[k] *= l_2k * std::exp(- l_2k * tau2[j]);
    }
    prob_j /= prob_j.sum();
    Rcpp::NumericVector prob_j_R = EigenVec_to_NumVec(prob_j);
    out[j] = int_sampling(K, 1, prob_j_R)[0];
  }
  
  return out;
}


/**
 * Function to perform "Step 2" in the blocked Gibbs update of the underlying  \ 
 * Dirichlet process mixture in the nonparametric Bayesian Lasso.              \
 *                                                                             \
 * Important: For computational efficiency, this function uses single          \ 
 * precision (floats)! If needed, the above function,                          \
 * "step2_blockedGibbs_dbl()", uses double precision.                          \
 *                                                                             \
 * @param tau2 VectorXd: Vector of size \eqn{p} with current values for tau2.  \
 * @param lambda2_all VectorXd: Vector of size \eqn{K} with all the current    \
 *  possible values for lambda2                                                \
 * @param omega VectorXd:VectorXd: Vector of size \eqn{K} with all the current \
 *  probabilities from each component.                                         \
 * @param K int: Truncation.                                                   \
 *                                                                             \
 * @return Rcpp::IntegerVector: Vector of updated cluster allocations          \
 */
// [[Rcpp::export]]
Rcpp::IntegerVector step2_blockedGibbs_flt(const Eigen::VectorXf & tau2,
                                           const Eigen::VectorXf & lambda2_all,
                                           const Eigen::VectorXf & omega,
                                           const int K) {
  
  // Pre-compute constants and prepare the returns
  const int p = tau2.size();
  Rcpp::IntegerVector out(p);
  
  for (int j = 0; j < p; ++j) {
    Eigen::VectorXf prob_j = omega; // Initialize vector of probabilities
    for (int k = 0; k < K; ++k) {
      float l_2k = lambda2_all[k] * 0.5;
      prob_j[k] *= l_2k * std::exp(- l_2k * tau2[j]);
    }
    prob_j /= prob_j.sum();
    Rcpp::NumericVector prob_j_R = EigenVecFloat_to_NumVec(prob_j);
    out[j] = int_sampling(K, 1, prob_j_R)[0];
  }
  
  return out;
}
