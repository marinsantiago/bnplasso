// -----------------------------------------------------------------------------
// Gibbs sampler that implements either:
//   (1) Nonparametric Bayesian Lasso.
//   (2) Bayesian Lasso.
//   (3) Bayesian Adaptive Lasso.
// -----------------------------------------------------------------------------

#include <RcppEigen.h>
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <chrono>
#include "utils.h"
#include "int_sampling.h"
#include "update_lambda_bnp.h"
#include "update_lambda_lasso.h"
#include "update_beta_lm.h"
#include "update_beta_spm.h"
#include "update_sigma2_lm.h"
#include "update_sigma2_spm.h"
#include "update_tau2.h"
#include "update_mu.h"
#include "mcmc_algorithm.h"


/**
 * Markov chain Monte Carlo (MCMC) algorithm to generate random draws from the \
 * posterior distribution of either: a nonparametric Bayesian Lasso, a         \ 
 * Bayesian Lasso, or a Bayesian Adaptive Lasso.                               \
 *                                                                             \
 * @return Rcpp::List: List containing the posterior draws as well as details  \
 * from the call of the function.                                              \
 */
// [[Rcpp::export]]
Rcpp::List mcmc_sampler(const Eigen::MatrixXd & X,
                        const Eigen::VectorXd & y,
                        const double a,
                        const double b,
                        const double alpha,
                        const bool intercept,
                        const std::string shrinkage_prior, 
                        const std::string variance_prior_type,
                        const bool polya_urn_sampler,
                        const bool sp_means,
                        const int max_iters,
                        const int burn_in,
                        const int thin,
                        const bool verbose) {
  
  // Pre-compute constants -----------------------------------------------------
  const Eigen::Index n = y.size();
  Eigen::MatrixXd Xt;
  Eigen::MatrixXd XtX;
  const double eps = 1e-06;
  Eigen::Index p;
  if (sp_means == true) {
    p = n;
  } else {
    p = X.cols();
    Xt = X.transpose();
    XtX = Xt * X; // Cross-product t(X) %*% X
    XtX = make_posdef(XtX);
  }
  
  // Initialize data structures to store posterior draws -----------------------
  int current_save = 0;
  // Number of draws after burn-in and thinning
  const int dd = std::floor((max_iters - burn_in) / thin);
  Eigen::VectorXd mu_out;
  if (intercept == true) { mu_out = Eigen::VectorXd::Zero(dd); }
  Eigen::VectorXd K_out(dd);
  Eigen::VectorXd sigma2_out(dd);
  Eigen::MatrixXd lambda2_out(dd, p);
  Eigen::MatrixXd beta_out(dd, p);
  Eigen::MatrixXd tau2_out(dd, p);
  Rcpp::IntegerMatrix clust_idx_out(dd, p);
  // Blocked Gibbs sampler parameters
  int K_trunc;
  Eigen::MatrixXd omega_out;
  if (polya_urn_sampler == false) {
    K_trunc = 250;
    omega_out = Eigen::MatrixXd::Zero(dd, K_trunc);
  }
  
  // Initialize model parameters -----------------------------------------------
  Rcpp::IntegerVector clust_indicators;
  Eigen::VectorXd lambda2(p);
  Eigen::VectorXd tau2(p);
  Eigen::VectorXd lambda2_candidates;
  Eigen::VectorXd lambda2_k;
  Eigen::VectorXd omega;
  Eigen::VectorXd beta;
  double sigma2 = 1.0;
  double mu = 0.0;
  int K_clust;
  if (shrinkage_prior == "bnplasso") {
    // Nonparametric Bayesian Lasso
    if (polya_urn_sampler == true) {
      // Generalized Polya urn sampling scheme
      const int n_clust_init = 4;
      // Randomly assign coeficients to some initial clusters
      clust_indicators = int_sampling(
        n_clust_init, p, Eigen::VectorXd::Ones(n_clust_init)
      );
      // Consider a sequence of lambda2_k values applying different shrinkage
      lambda2_k = Eigen::VectorXd::Zero(n_clust_init);
      lambda2_k[0] = 0.01;
      lambda2_k[1] = 0.1;
      lambda2_k[2] = 0.5;
      lambda2_k[3] = 200;
      for (Eigen::Index jj = 0; jj < p; ++jj) {
        // Get the corresponding exponential mixing rates
        const double mix_rates = lambda2_k[clust_indicators[jj] - 1];
        // Draw each tau_j from their respective prior distributions
        tau2[jj] = std::max(R::rexp(1.0 / (mix_rates / 2.0)), eps);
      }
    } else {
      // Blocked Gibbs sampler
      lambda2_candidates = Eigen::VectorXd::Zero(K_trunc);
      omega = Eigen::VectorXd::Zero(K_trunc);
      Eigen::VectorXd nu(K_trunc - 1);
      double remaining = 1.0; // prod_{l<k} (1 - nu_l)
      for (int kk = 0; kk < (K_trunc - 1); ++kk) {
        // Sample lambda2 candidates from G_0
        lambda2_candidates[kk] = std::max(R::rgamma(a, 1.0 / b), eps);
        // Sample nu
        const double nu = R::rbeta(1.0, alpha);
        // Sample omega
        omega[kk] = nu * remaining;
        remaining *= (1.0 - nu);
      }
      lambda2_candidates[K_trunc - 1] = std::max(R::rgamma(a, 1.0 / b), eps);
      omega[K_trunc - 1] = remaining;
      // Sample the cluster allocations given the above omega
      clust_indicators = int_sampling(K_trunc, p, omega);
      // Sample tau2
      for (Eigen::Index jj = 0; jj < p; ++jj) {
        // Get the corresponding exponential mixing rates
        const double mix_rates = lambda2_candidates[clust_indicators[jj] - 1];
        // Draw each tau_j from their respective prior distributions
        tau2[jj] = std::max(R::rexp(1.0 / (mix_rates / 2.0)), eps);
      }
    }
  } else if (shrinkage_prior == "b_lasso") {
    // Bayesian Lasso: Only one lambda
    const double lambda2_single = std::max(R::rgamma(a, 1.0 / b), eps);
    for (Eigen::Index jj = 0; jj < p; ++jj) {
      lambda2[jj] = lambda2_single;
      // Draw tau2 from its respective prior distribution
      tau2[jj] = std::max(R::rexp(1.0 / (lambda2_single / 2.0)), eps);
    }
  } else if (shrinkage_prior == "b_adapt_lasso") {
    // Bayesian adaptive Lasso: p lambdas
    for (Eigen::Index jj = 0; jj < p; ++jj) {
      const double lambda2_j = std::max(R::rgamma(a, 1.0 / b), eps);
      lambda2[jj] = lambda2_j;
      // Draw tau2 from its respective prior distribution
      tau2[jj] = std::max(R::rexp(1.0 / (lambda2_j / 2.0)), eps);
    }
  }
  
  // MCMC algorithm ------------------------------------------------------------
  const auto start = std::chrono::steady_clock::now(); // Start wall-clock time
  Progress progress(max_iters, verbose);
  for (int iter = 0; iter < max_iters; ++iter) {
    if (Progress::check_abort()) { Rcpp::stop("Algorithm halted."); }
    // Update lambda2
    if (shrinkage_prior == "bnplasso") {
      // Nonparametric Bayesian Lasso 
      if (polya_urn_sampler == true) {
        // Generalized Polya urn sampling scheme
        const Rcpp::List out_DP_Polya = update_lambda_DP_Polya_urn(
          tau2, lambda2_k, clust_indicators, a, b, alpha
        );
        clust_indicators = out_DP_Polya["clust_indicators_out"];
        lambda2_k = out_DP_Polya["lambda2_k_out"];
        K_clust = out_DP_Polya["K"];
        for (Eigen::Index jj = 0; jj < p; ++jj) {
          lambda2[jj] = lambda2_k[clust_indicators[jj] - 1];
        }
      } else {
        // Blocked Gibbs sampler 
        const Rcpp::List out_DP_blocked_Gibbs = update_lambda_DP_blocked_Gibbs(
          tau2, lambda2_candidates, clust_indicators, omega, a, b, alpha
        );
        lambda2_candidates = out_DP_blocked_Gibbs["lambda2_candidates_out"];
        clust_indicators = out_DP_blocked_Gibbs["clust_indicators_out"];
        omega = out_DP_blocked_Gibbs["omega_out"];
        K_clust = out_DP_blocked_Gibbs["K"];
        for (Eigen::Index jj = 0; jj < p; ++jj) {
          lambda2[jj] = lambda2_candidates[clust_indicators[jj] - 1];
        }
      }
    } else if (shrinkage_prior == "b_lasso") {
      // Bayesian Lasso
      lambda2 = update_lambda2_Bayesian_lasso(tau2, a, b);
    } else if (shrinkage_prior == "b_adapt_lasso") {
      // Bayesian adaptive Lasso
      lambda2 = update_lambda2_Bayesian_adapt_lasso(tau2, a, b);
    }
    // Update beta, sigma2, and tau2
    if (variance_prior_type == "independent") {
      // Assuming an independent variance prior
      if (sp_means == false) {
        // Linear regression problem
        beta = update_beta_lm_ind_sigma(Xt, XtX, y, tau2, sigma2, mu);
        sigma2 = update_sigma2_ind_prior(X, y, beta, mu, intercept);
      } else {
        // Sparse means problem
        beta = update_beta_spm_ind_sigma(y, tau2, sigma2);
        sigma2 = update_sigma2_spm_ind_prior(y, beta);
      }
      tau2 = update_tau2_ind_sigma(beta, lambda2);
    } else if (variance_prior_type == "conjugate") {
      // Assuming a conjugate variance prior
      if (sp_means == false) {
        // Linear regression problem
        beta = update_beta_lm_conj_sigma(Xt, XtX, y, tau2, sigma2, mu);
        sigma2 = update_sigma2_conj_prior(X, y, beta, tau2, mu, intercept);
      } else {
        // Sparse means problem
        beta = update_beta_spm_conj_sigma(y, tau2, sigma2);
        sigma2 = update_sigma2_spm_conj_prior(y, beta, tau2);
      }
      tau2 = update_tau2_conj_sigma(beta, lambda2, sigma2);
    }
    // Update mu
    if (intercept == true) { mu = update_mu(X, y, beta, sigma2); }
    // Store draws after burn-in and thinning
    if (iter >= burn_in) {
      if (((iter - burn_in + 1) % thin == 0)) {
        sigma2_out[current_save] = sigma2;
        lambda2_out.row(current_save) = lambda2;
        beta_out.row(current_save) = beta;
        tau2_out.row(current_save) = tau2;
        if (intercept == true) { mu_out[current_save] = mu; }
        if (shrinkage_prior == "bnplasso") {
          K_out[current_save] = K_clust;
          clust_idx_out.row(current_save) = clust_indicators;
          if (polya_urn_sampler == false) {
            omega_out.row(current_save) = omega;
          }
        }
        ++current_save;
      }
    }
    progress.increment();
  }
  const auto end = std::chrono::steady_clock::now(); // End wall-clock time
  const double elapsed = std::chrono::duration<double>(end - start).count();
  
  // Prepare the returns -------------------------------------------------------
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("post.lambda2") = lambda2_out,
    Rcpp::Named("post.sigma2") = sigma2_out,
    Rcpp::Named("post.beta")  = beta_out,
    Rcpp::Named("post.tau2") = tau2_out,
    Rcpp::Named("elapsed") = elapsed
  );
  if (intercept == true) { out["post.mu"] = mu_out; }
  if (shrinkage_prior == "bnplasso") {
    out["post.K"] = K_out;
    out["post.clust_idx"] = clust_idx_out;
    if (polya_urn_sampler == false) { out["post.omega"] = omega_out; }
  }
  
  return out;
}
