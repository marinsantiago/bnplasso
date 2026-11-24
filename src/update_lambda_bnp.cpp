// -----------------------------------------------------------------------------
// Update lambda under a Dirichlet process prior.
// This can be done through either:
//   (a) A generalized Polya urn scheme.
//   (b) A blocked Gibbs sampler.
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "utils.h"
#include "int_sampling.h"
#include "update_lambda_bnp.h"


/**
 * Update the lambda parameters under a Dirichlet process prior.               \
 * This function employs a generalized Polya urn scheme.                       \
 *                                                                             \
 * @param tau2 Eigen::VectorXd: Current value of the tau2 parameters.          \
 * @param lambda2_k Eigen::VectorXd: Current value of lambda2_k.               \
 * @param clust_indicators Rcpp::IntegerVector: Current cluster labels.        \
 * @param a double: Shape parameter of the prior centering measure.            \
 * @param b double: Rate parameter of the prior centering measure.             \
 * @param alpha double: Concentration parameter of the DP prior.               \
 *                                                                             \
 * @return Rcpp::List: Updated values of lambda2_k and the cluster labels.     \
 */
// [[Rcpp::export]]
Rcpp::List update_lambda_DP_Polya_urn(
    const Eigen::VectorXd & tau2,
    const Eigen::VectorXd & lambda2_k,
    const Rcpp::IntegerVector & clust_indicators,
    const double a, 
    const double b,
    const double alpha
) {
  
  // Pre-compute constants -----------------------------------------------------
  const Eigen::Index p = tau2.size();
  const int K_in = lambda2_k.size();
  const double eps = 1e-08;
  
  // Compact the clsuter labels and transform them into 0-based indices --------
  std::vector<double> lambda2;
  std::vector<int> new_label(K_in, -1);
  std::vector<int> z(p);
  std::vector<int> n;
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    const int current_clust = clust_indicators[jj] - 1; // C++ indexing 
    if (new_label[current_clust] < 0) {
      new_label[current_clust] = static_cast<int>(lambda2.size());
      lambda2.push_back(lambda2_k[current_clust]);
      n.push_back(0);
    }
    z[jj] = new_label[current_clust];
    ++n[z[jj]];
  }
  int K = static_cast<int>(lambda2.size());
  
  // Update cluster membership indicators using a CRP --------------------------
  const double prob_new_const = 0.5 * alpha * a * std::pow(b, a);
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    const double current_tau2 = tau2[jj];
    const int current_label = z[jj];
    // Remove tau_j; drop its cluster if now empty (move last cluster into slot)
    if (--n[current_label] == 0) {
      const int last_idx = K - 1; // C++ indexing 
      if (current_label != last_idx) {
        lambda2[current_label] = lambda2[last_idx];
        n[current_label] = n[last_idx];
        for (Eigen::Index ii = 0; ii < p; ++ii) {
          if (z[ii] == last_idx) { z[ii] = current_label; }
        }
      }
      lambda2.pop_back();
      n.pop_back();
      --K;
    }
    // Probabilities of joining an existing cluster or joining a new one
    Eigen::VectorXd probs(K + 1);
    // Existing clusters
    for (int kk = 0; kk < K; ++kk) {
      const double lambda2_half = 0.5 * lambda2[kk];
      probs[kk] = static_cast<double>(n[kk]) * lambda2_half;
      probs[kk] *= std::exp(-lambda2_half * current_tau2);
    }
    // New cluster 
    probs[K] = prob_new_const * std::pow(b + current_tau2 / 2, -(a + 1));
    // Sample the updated cluster label ----------------------------------------
    const int label_new = int_sampling(K + 1, 1, probs)[0] - 1; // C++ indexing 
    if (label_new == K) {
      // New cluster: Sample lambda2 from its posterior given tau_j only
      lambda2.push_back(
        std::max(R::rgamma(a + 1.0, 1.0 / (b + 0.5 * current_tau2)), eps)
      );
      n.push_back(1);
      ++K;
    } else {
      ++n[label_new];
    }
    z[jj] = label_new;
  }
  
  // Update the lambda2's from their full conditional distributions  -----------
  Eigen::VectorXd sums_tau = Eigen::VectorXd::Zero(K);
  for (Eigen::Index jj = 0; jj < p; ++jj) { sums_tau[z[jj]] += tau2[jj]; }
  for (int kk = 0; kk < K; ++kk) {
    lambda2[kk] = std::max(
      R::rgamma(a + n[kk], 1.0 / (b + 0.5 * sums_tau[kk])), eps
    );
  }
  
  // Prepare the returns -------------------------------------------------------
  Rcpp::IntegerVector clust_indicators_out(p);
  Eigen::VectorXd lambda2_k_out(K); 
  for (int kk = 0; kk < K; ++kk) { lambda2_k_out[kk] = lambda2[kk]; }
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    clust_indicators_out[jj] = z[jj] + 1; // One-based indexing 
  }
  
  return Rcpp::List::create(
    Rcpp::Named("clust_indicators_out") = clust_indicators_out,
    Rcpp::Named("lambda2_k_out") = lambda2_k_out,
    Rcpp::Named("K") = K
  );
}


/**
 * Update the lambda parameters under a Dirichlet process prior.               \
 * This function employs a blocked Gibbs sampler.                              \
 *                                                                             \
 * @param tau2 Eigen::VectorXd: Current value of the tau2 parameters.          \
 * @param lambda2_candidates Eigen::VectorXd: Current candidate values of      \
 *    lambda2_k.                                                               \      
 * @param clust_indicators Rcpp::IntegerVector: Current cluster labels.        \
 * @param omega Eigen::VectorXd: Current value of the stick-breaking weights.  \
 * @param a double: Shape parameter of the prior centering measure.            \
 * @param b double: Rate parameter of the prior centering measure.             \
 * @param alpha double: Concentration parameter of the DP prior.               \
 *                                                                             \
 * @return Rcpp::List: Updated values of lambda2_k and the cluster labels.     \
 */
// [[Rcpp::export]]
Rcpp::List update_lambda_DP_blocked_Gibbs(
    const Eigen::VectorXd & tau2,
    const Eigen::VectorXd & lambda2_candidates,   
    const Rcpp::IntegerVector & clust_indicators, 
    const Eigen::VectorXd & omega,
    const double a,
    const double b,
    const double alpha
) {
  
  // Pre-compute constants -----------------------------------------------------
  const Eigen::Index p = tau2.size();
  const int K = lambda2_candidates.size();
  const double eps = 1e-08;
  
  // Cluster sizes and tau2 sums under the current allocations -----------------
  std::vector<int> n(K, 0);
  std::vector<double> sum_tau2(K, 0.0);
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    const int k = clust_indicators[jj] - 1; // C++ indexing 
    ++n[k];
    sum_tau2[k] += tau2[jj];
  }
  
  // Update lambda2_candidates -------------------------------------------------
  Eigen::VectorXd lambda2_candidates_out(K);
  for (int kk = 0; kk < K; ++kk) {
    // Gamma(a + n_k, b + sum/2); reduces to the prior G0 when n_k = 0.
    lambda2_candidates_out[kk] = std::max(
      R::rgamma(a + n[kk], 1.0 / (b + 0.5 * sum_tau2[kk])), eps
    );
  }
  
  // Update cluster allocations ------------------------------------------------
  Rcpp::IntegerVector clust_indicators_out(p);
  Eigen::VectorXd lambda2_half = 0.5 * lambda2_candidates_out;
  // Compute omega_k * lambda2_k / 2
  Eigen::VectorXd coef = omega.cwiseProduct(lambda2_half); 
  Eigen::VectorXd probs(K);
  std::fill(n.begin(), n.end(), 0); // Reuse to update the new counts
  for (Eigen::Index jj = 0; jj < p; ++jj) {
    const double t2 = tau2[jj];
    for (int kk = 0; kk < K; ++kk) {
      probs[kk] = coef[kk] * std::exp(-lambda2_half[kk] * t2);
    }
    // Sample the cluster label
    const int s = int_sampling(K, 1, probs)[0];
    clust_indicators_out[jj] = s;
    ++n[s - 1]; // C++ indexing
  }
  
  // Update the stick-breaking weights -----------------------------------------
  Eigen::VectorXd omega_out(K);
  int n_above = static_cast<int>(p); // number of tau's with label > k
  double remaining = 1.0; // prod_{l<k} (1 - nu_l)
  for (int kk = 0; kk < (K - 1); ++kk) {
    n_above -= n[kk];
    // Sample nu
    const double nu = R::rbeta(1.0 + n[kk], alpha + n_above);
    omega_out[kk] = nu * remaining;
    remaining *= (1.0 - nu);
  }
  omega_out[K - 1] = remaining;
  
  // Prepare the returns -------------------------------------------------------
  return Rcpp::List::create(
    Rcpp::Named("lambda2_candidates_out") = lambda2_candidates_out,
    Rcpp::Named("clust_indicators_out") = clust_indicators_out,
    Rcpp::Named("omega_out") = omega_out,
    Rcpp::Named("K") = K
  );
}
