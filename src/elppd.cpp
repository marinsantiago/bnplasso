// -----------------------------------------------------------------------------
// Expected log pointwise predictive density (elppd)
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "elppd.h"


//' Expected log pointwise predictive density (elppd)
//'
//' This function computes the expected log pointwise predictive density (elppd) 
//' for an object of class \code{lmBayes}. Given a test (held-out) set, 
//' \eqn{\{y_{\text{test},i},\,
//' \mathbf{x}_{\text{test},i}\}_{i=1}^{n_{\text{test}}}}, the elppd is defined 
//' as 
//' \deqn{\text{elppd} = \frac{1}{n_{\text{test}}}\sum_{i=1}^{n_{\text{test}}}
//' \log\left(\frac{1}{S}\sum_{s=1}^{S}p\left({y}_{\text{test}, i} | 
//' \mathbf{x}_{\text{test}, i}^{_{'}}\boldsymbol{\beta}_{(s)}, 
//' \sigma^{2}_{(s)}\right)\right),}
//' where \eqn{p()} is the likelihood of the observed data and 
//' \eqn{\boldsymbol{\beta}_{(s)}} and \eqn{\sigma^{2}_{(s)}} are the 
//' \eqn{s}-th draws of the coefficient vector and the sampling variance in the 
//' MCMC algorithm, respectively. 
//'
//' @param object An object of class \code{lmBayes}.
//' @param X_new A matrix of predictors from the test data, where each row is an 
//'   observation. 
//' @param y_new A vector of responses from the test data of size 
//'   \eqn{n_{\text{test}}}.
//' 
//' @return The elppd as a numeric scalar.
//'
//' @references
//'
//' A. Gelman, J. Carlin, H. Stern, D. Dunson, and A. Vehtari (2013), 
//' \emph{Bayesian Data Analysis}, 3. Chapman & Hall/CRC.
//' 
//' @author Santiago Marin
//'
// [[Rcpp::export]]
double elppd(const Rcpp::List & object, 
             const Eigen::MatrixXd & X_new,
             const Eigen::VectorXd & y_new) {
  
  // Input validation ----------------------------------------------------------
  const bool is_lmBayes = object.inherits("lmBayes");
  const int n_draws = object["n.draws"];
  const int n_preds = object["n.preds"];
  if (y_new.array().hasNaN() == true) {
    Rcpp::stop("y_new is not a valid response vector");
  }
  const int n_test = y_new.size();
  if (X_new.array().hasNaN() == true) {
    Rcpp::stop("X_new is not a valid matrix of covariates");
  }
  if (n_test != X_new.rows()) {
    Rcpp::stop("Incorrect number of rows in X_new");
  }
  if (n_preds != X_new.cols()) {
    Rcpp::stop("Incorrect number of columns in X_new");
  }

  // Compute elppd -------------------------------------------------------------
  static const double log2pi = std::log(2.0 * M_PI);
  const double eps = 1e-16;
  const Eigen::MatrixXd post_beta = object["post.beta"];
  // Pre-compute all linear predictors (without the intercept) for the new data
  Eigen::MatrixXd lin_preds = X_new * post_beta.transpose();
  // Add the intercept if needed
  const bool intercept = object["intercept"];
  if (intercept == true) {
    const Eigen::VectorXd mu = object["post.mu"];
    lin_preds.rowwise() += mu.transpose();
  }
  Eigen::VectorXd post_sigma2 = object["post.sigma2"];
  post_sigma2.array() += eps;
  // Log point-wise predictive density
  Eigen::VectorXd lppd(n_test);
  for (int ii = 0; ii < n_test; ++ii) {
    // Compute the normal density
    const Eigen::VectorXd lin_pred_i = lin_preds.row(ii);
    Eigen::VectorXd exponent = y_new[ii] - lin_pred_i.array();
    exponent = exponent.array().square();
    exponent = -0.5 * (exponent.array() / post_sigma2.array());
    Eigen::VectorXd denom = post_sigma2.array().sqrt().log(); 
    denom.array() += 0.5 * log2pi;
    Eigen::VectorXd dens = exponent - denom;
    dens = dens.array().exp() + 1e-25;
    double lppd_i = dens.mean();
    lppd_i = std::log(lppd_i);
    lppd[ii] = lppd_i;
  }
  return lppd.mean();
}
