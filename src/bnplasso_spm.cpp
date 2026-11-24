// -----------------------------------------------------------------------------
// Nonparametric Bayesian Lasso for sparse means problems
// -----------------------------------------------------------------------------
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "utils.h"
#include "mcmc_algorithm.h"
#include "bnplasso_spm.h"


//' Fit sparse means models with a nonparametric Bayesian Lasso prior
//'
//' This function fits sparse means models with a nonparametric Bayesian 
//' Lasso prior (Marin et al., 2026).
//'
//' @param y Response variable. A numeric vector size \eqn{n}.
//' @param shrinkage_prior A character string denoting which type of shrinkage 
//'   prior should be employed. The options are either: (1) \code{"bnplasso"} 
//'   which implements the nonparametric Bayesian Lasso as in 
//'   Marin et al. (2026), (2) \code{"b_lasso"} which implements the Bayesian 
//'   Lasso as in Park and Casella (2008), or (3) \code{"b_adapt_lasso"} which 
//'   implements the Bayesian adaptive Lasso as in Leng et al. (2014). 
//'   Default is \code{"bnplasso"}.
//' @param a A positive scalar. If \code{prior = "bnplasso"}, it corresponds to
//'   the \bold{shape} parameter in the gamma distribution used as a centering 
//'   measure in the DP prior. If \code{prior = "b_lasso"}, it corresponds to 
//'   the \bold{shape} parameter in the gamma distribution used as prior for the 
//'   shrinkage parameter, \eqn{\lambda^{2}}. If \code{prior = "b_adapt_lasso"}, 
//'   it corresponds to the \bold{shape} parameter in the gamma distribution 
//'   used as prior for the shrinkage parameters, \eqn{\lambda_{j}^{2}}, for 
//'   \eqn{j\in\{1,\dots,p\}}. If not provided, the function will attempt to 
//'   determine an appropriate value.
//' @param b A positive scalar. If \code{prior = "bnplasso"}, it corresponds to 
//'   the \bold{rate} parameter in the gamma distribution used as a centering 
//'   measure in the DP prior. If \code{prior = "b_lasso"}, it corresponds to 
//'   the \bold{rate} parameter in the gamma distribution used as prior for the 
//'   shrinkage parameter, \eqn{\lambda^{2}}. If \code{prior = "b_adapt_lasso"}, 
//'   it corresponds to the \bold{rate} parameter in the gamma distribution 
//'   used as prior for the shrinkage parameters, \eqn{\lambda_{j}^{2}}, for 
//'   \eqn{j\in\{1,\dots,p\}}. If not provided, the function will attempt to 
//'   determine an appropriate value.
//' @param alpha A positive scalar. If \code{prior = "bnplasso"}, it 
//'   corresponds to the \bold{concentration} parameter in the DP prior. If 
//'   \code{prior = "b_lasso"} or \code{prior = "b_adapt_lasso"},
//'   the argument is ignored. If not provided, the function will attempt to 
//'   determine an appropriate value.
//' @param variance_prior_type A character string denoting whether the prior 
//'   on the sampling variance should be an independent-type prior or a 
//'   conjugate-type prior. See Moran et al. (2019) for details. The options are 
//'   either \code{"independent"} or \code{"conjugate"}. If not provided, the 
//'   function will attempt to determine an appropriate prior.
//' @param max_iters A positive integer corresponding to the total number of 
//'   MCMC iterations. Default is 10000.
//' @param burn_in A positive integer corresponding to the number of draws 
//'   discarded as burn-in. It should be smaller than \code{max_iters}.
//'   Default is 20000.
//' @param thin A positive integer specifying the period for saving samples.
//'   Default is 2.
//' @param polya_urn_sampler logical. Whether a generalized Pólya urn sampling 
//'   scheme or a blocked Gibbs sampling scheme should be employed to update 
//'   the Dirichlet process mixture in the MCMC algorithm. If \code{TRUE}, a 
//'   generalized Pólya urn sampling scheme would be used. Default is 
//'   \code{TRUE}. This argument is ignored if \code{prior != "bnp.lasso"}
//' @param verbose Logical. Whether a progress bar should be printed in the
//'   the console. Default is \code{TRUE}.
//'
//' @return An object of S3 class, \code{'spmBayes'}, containing:
//' \itemize{
//'   \item \code{post.beta}: A matrix of size \code{n.draws}-by-\code{n.preds}, 
//'   where each row is a posterior draw of the regression coefficients.
//'   \item \code{post.sigma2}: A vector of size \code{n.draws}, where each 
//'   element is a posterior draw of the sampling variance. 
//'   \item \code{post.tau2}: A matrix of size \code{n.draws}-by-\code{n.preds}, 
//'   where each row is a posterior draw of the latent parameters 
//'   \eqn{\tau_{j}^{2}}.
//'   \item \code{post.lambda2}: A matrix of size 
//'   \code{n.draws}-by-\code{n.preds}, where each row is a posterior draw
//'   of the shrinkage parameters \eqn{\lambda_{j}^{2}}.
//'   \item \code{post.clust_idx}: If \code{prior = "bnp.lasso"}, a matrix of 
//'   size \code{n.draws}-by-\code{n.preds}, where each row indicates to which 
//'   cluster the regression coefficients belong to.
//'   \item \code{post.K}: If \code{prior = "bnp.lasso"}, a vector of size 
//'   \code{n.draws}, where each element indicates the number of clusters in the 
//'   corresponding MCMC iteration.
//'   \item \code{elapsed}: The elapsed (wall-clock) time of the MCMC sampler, 
//'   in seconds.  
//'   \item \code{a}: If \code{prior = "bnp_lasso"}, the \bold{shape} parameter 
//'   in the gamma distribution used as a centering measure in the DP prior. 
//'   If \code{prior = "b_lasso"}, the \bold{shape} parameter in the gamma 
//'   distribution used as prior in the shrinkage parameter, \eqn{\lambda^{2}}. 
//'   If \code{prior = "b_adapt_lasso"}, the \bold{shape} parameter in the gamma 
//'   distribution used as prior in the shrinkage parameters, 
//'   \eqn{\lambda_{j}^{2}}, for \eqn{j\in\{1,\dots,p\}}.
//'   \item \code{b}: If \code{prior = "bnplasso"}, the \bold{rate} parameter 
//'   in the gamma distribution used as a centering measure in the DP prior. 
//'   If \code{prior = "b_lasso"}, the \bold{rate} parameter in the gamma 
//'   distribution used as prior in the shrinkage parameter, \eqn{\lambda^{2}}. 
//'   If \code{prior = "b_adapt_lasso"}, the \bold{rate} parameter in the gamma 
//'   distribution used as prior in the shrinkage parameters, 
//'   \eqn{\lambda_{j}^{2}}, for \eqn{j\in\{1,\dots,p\}}.
//'   \item \code{alpha}: If \code{prior = "bnplasso"}, the \bold{concentration} 
//'   parameter in the Dirichlet process prior.
//'   \item \code{variance.prior.type}: Whether the variance on the sampling 
//'   variance was an independent-type prior or a conjugate-type prior.
//'   \item \code{max.iters}: The total number of MCMC iterations.
//'   \item \code{burn.in}: The number of draws discarded as burn-in.
//'   \item \code{thin}: The period for saving draws.
//'   \item \code{n.obs}: The sample size.   
//'   \item \code{n.preds}: The number of predictors.
//'   \item \code{n.draws}: The number of posterior draws after burn-in and
//'   thinning.
//'   \item \code{y}: Vector of responses.
//'   \item \code{loglik}: Matrix of size \code{n.draws}-by-\code{n.obs} with 
//'   the log-likelihood of each observation at each MCMC iteration.
//'   \item \code{post.pred.fitted.values}: A matrix of size 
//'   \code{n.draws}-by-\code{n.obs}, where each row is a draw from the 
//'   posterior predictive distribution of the fitted values.
//'   \item \code{post.pred.residuals}: A matrix of size 
//'   \code{n.draws}-by-\code{n.obs}, where each row is a draw from the 
//'   posterior predictive distribution of the residuals.  
//' }
//' 
//' @references
//'
//' C. Leng, MN. Tran, and D. Nott (2014), Bayesian adaptive Lasso. 
//' \emph{Ann Inst Stat Math}, 66:221-244
//'
//' S. Marin, B. Long,and A. H. Westveld (2026), Adaptive Shrinkage with a 
//' Nonparametric Bayesian Lasso.\emph{Journal of Computational and Graphical 
//' Statistics}, 35(2):854-864. doi:10.1080/10618600.2025.2572327
//' 
//' G. E. Moran, V. Rockova, and E. I. George (2019), Variance Prior Forms for 
//' High-Dimensional Bayesian Variable Selection. \emph{Bayesian Analysis},
//' 14(4):1091-1119.
//' 
//' T. Park, and G. Casella (2008), The Bayesian Lasso. 
//' \emph{Journal of the American Statistical Association}, 103(482):681-686.
//' 
//' @author Santiago Marin
//'
//' @encoding UTF-8
//'
// [[Rcpp::export]]
Rcpp::List bnplasso_spm(
    const Eigen::VectorXd & y,
    Rcpp::Nullable<Rcpp::String> shrinkage_prior = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> a = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> b = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> alpha = R_NilValue,
    Rcpp::Nullable<Rcpp::String> variance_prior_type = R_NilValue,
    const int max_iters = 10000,
    const int burn_in = 2000,
    const int thin = 2,
    const bool polya_urn_sampler = true,
    const bool verbose = true
) {
  
  // Input validation ----------------------------------------------------------
  if (y.array().hasNaN() == true) {
    Rcpp::stop("y is not a valid response vector");
  }
  const int n = y.size();
  if (shrinkage_prior.isNotNull()) {
    const std::string s = Rcpp::as<std::string>(shrinkage_prior);
    if (s != "bnplasso" && s != "b_lasso" && s != "b_adapt_lasso") {
      Rcpp::stop("Incorrect shrinkage prior");
    }
  }
  if (variance_prior_type.isNotNull()) {
    std::string v = Rcpp::as<std::string>(variance_prior_type);
    if (v != "independent" && v != "conjugate") {
      Rcpp::stop("Incorrect variance prior type");
    }
  }
  if (a.isNotNull()) {
    double aa = Rcpp::as<double>(a);
    if (!(aa > 0)) { Rcpp::stop("a is not a postive scalar"); }
  }
  if (b.isNotNull()) {
    double bb = Rcpp::as<double>(b);
    if (!(bb > 0)) { Rcpp::stop("b is not a postive scalar"); }
  }
  if (alpha.isNotNull()) {
    double aalpha = Rcpp::as<double>(alpha);
    if (!(aalpha > 0)) { Rcpp::stop("alpha is not a postive scalar"); }
  }
  if (max_iters <= 0) { Rcpp::stop("max_iters is not a postive integer"); }
  if (burn_in <= 0) { Rcpp::stop("burn_in is not a postive integer"); }
  if (thin <= 0) { Rcpp::stop("thin is not a postive integer"); }
  if (burn_in >= max_iters) {
    Rcpp::stop("burn_in should be smaller than max_iters"); 
  }
  if (thin >= max_iters) { 
    Rcpp::stop("thin should be smaller than max_iters");
  }
  
  // Set default hyper-parameters (if not provided) ----------------------------
  const bool sp_means = true; // Sparse means problems
  const bool intercept = false; 
  const Eigen::MatrixXd X = Eigen::MatrixXd::Zero(1, 1);
  std::string sshrinkage_prior = "bnplasso"; // Default
  if (shrinkage_prior.isNotNull()) {
    sshrinkage_prior = Rcpp::as<std::string>(shrinkage_prior);
  }
  double aa = 0.1; // Default
  if (a.isNotNull()) { aa = Rcpp::as<double>(a); }
  double bb = std::max(0.01 / n, 1e-10); // Default
  if (sshrinkage_prior != "bnplasso") { bb = 0.01; }
  if (b.isNotNull()) { bb = Rcpp::as<double>(b); }
  // Choose alpha so that prior K is approx. 3
  const int p = y.size();
  double aalpha = find_alpha(p, 3.0);
  if (sshrinkage_prior != "bnplasso") { aalpha = 1.0; }
  if (alpha.isNotNull()) { aalpha = Rcpp::as<double>(alpha); }
  std::string vvariance_prior_type = "independent"; // Default
  if (variance_prior_type.isNotNull()) {
    vvariance_prior_type = Rcpp::as<std::string>(variance_prior_type);
  }
  
  // Run Gibbs sampler ---------------------------------------------------------
  Rcpp::List bnplasso_out = mcmc_sampler(
    X, y, aa, bb, aalpha, intercept, sshrinkage_prior, vvariance_prior_type,
    polya_urn_sampler, sp_means, max_iters, burn_in, thin, verbose
  );
  
  // Prepare the returns -------------------------------------------------------
  const int n_draws = std::floor((max_iters - burn_in) / thin);
  bnplasso_out["n.obs"] = n;
  bnplasso_out["n.preds"] = p;
  bnplasso_out["n.draws"] = n_draws;
  bnplasso_out["a"] = aa;
  bnplasso_out["b"] = bb;
  if (sshrinkage_prior == "bnplasso") { bnplasso_out["alpha"] = aalpha; }
  bnplasso_out["shrinakge.prior"] = sshrinkage_prior;
  bnplasso_out["variance.prior.type"] = vvariance_prior_type;
  bnplasso_out["max.iters"] = max_iters;
  bnplasso_out["burn.in"] = burn_in;
  bnplasso_out["thin"] = thin;
  
  // Posterior predictive fitted values and residuals --------------------------
  const Eigen::MatrixXd post_beta = bnplasso_out["post.beta"];
  Eigen::VectorXd post_sigma = bnplasso_out["post.sigma2"];
  post_sigma.array() += 1e-16;
  post_sigma = post_sigma.array().sqrt();
  // Note: In "post_pred_fits" and "post_pred_res", each row corresponds to an
  // MCMC draw and each column to an observation.
  Eigen::MatrixXd post_pred_fits(n_draws, n);
  Eigen::MatrixXd post_pred_res(n_draws, n);
  Eigen::MatrixXd loglik(n_draws, n);
  for (int ss = 0; ss < n_draws; ++ss) {
    const Eigen::VectorXd beta_post_s = post_beta.row(ss);
    const double sigma_s = post_sigma[ss];
    Eigen::VectorXd epsilon(n); // Gaussian random noise
    Eigen::VectorXd log_dens(n);
    for (int ii = 0; ii < n; ++ii) {
      epsilon[ii] = R::rnorm(0.0, sigma_s);
      log_dens[ii] = R::dnorm4(y[ii], beta_post_s[ii], sigma_s, true);
    }
    const Eigen::VectorXd fit_val_s = beta_post_s + epsilon;
    post_pred_fits.row(ss) = fit_val_s;
    post_pred_res.row(ss) = y - fit_val_s;
    loglik.row(ss) = log_dens;
  }
  bnplasso_out["y"] = y;
  bnplasso_out["loglik"] = loglik;
  bnplasso_out["post.pred.fitted.values"] = post_pred_fits;
  bnplasso_out["post.pred.residuals"] = post_pred_res;
  bnplasso_out.attr("class") = "spmBayes";
  return bnplasso_out;
}
