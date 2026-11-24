// -----------------------------------------------------------------------------
// Various helper routines
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include "utils.h"


/**
 * Identify whether or not the integer x is within the integer vector v.       \ 
 *                                                                             \
 * @param x int: Integer value.                                                \
 * @param v std::vector<int>: Integer vector.                                  \
 *                                                                             \
 * @return bool: Whether the integer x is within the integer vector v.         \
 */
// [[Rcpp::export]]
bool within_integer_vec(const int x, const std::vector<int> v) {
  const int n = v.size();
  for (int ii = 0; ii < n; ++ii) {
    if (v[ii] == x) { return true; }
  }
  return false;
}


/**
 * Centers the columns of a numeric matrix.                                    \
 *                                                                             \
 * @param x Eigen::MatrixXd. Matrix to be centered.                            \
 *                                                                             \
 * @return Eigen::MatrixXd: Centered matrix.                                   \
 */
// [[Rcpp::export]]
Eigen::MatrixXd center(const Eigen::MatrixXd & x) { 
  Eigen::MatrixXd x_out = x;
  // Column means
  const Eigen::RowVectorXd means = x.colwise().mean();
  // Center each column
  x_out.rowwise() -= means;
  return x_out;
}


/**
 * Make a matrix m positive definite.                                          \
 */
// [[Rcpp::export]]
Eigen::MatrixXd make_posdef(const Eigen::MatrixXd & m) {
  const Eigen::Index d = m.rows(); // dim
  if (d != m.cols()) { Rcpp::stop("m must be a square matrix."); }
  if (d == 0) { return Eigen::MatrixXd(0, 0); }
  if (!m.allFinite()) { Rcpp::stop("m must contain only finite values."); }
  const Eigen::MatrixXd m_sym = 0.5 * (m + m.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(m_sym);
  if (eig.info() != Eigen::Success) { Rcpp::stop("Eigendecomposition failed"); }
  Eigen::VectorXd evals = eig.eigenvalues();
  const double eps = std::numeric_limits<double>::epsilon();
  const double dd = static_cast<double>(d);
  const double tol = 2.0 * dd * evals.cwiseAbs().maxCoeff() * eps; 
  evals = evals.cwiseMax(tol);
  const Eigen::MatrixXd & evecs = eig.eigenvectors();
  return evecs * evals.asDiagonal() * evecs.transpose();
}


/**
 * Given a number of covariates, this functions returns the concentration      \
 *  parameter in a DP (alpha), so that the prior number of clusters is the     \
 *  closest to "target".                                                       \
 */
// [[Rcpp::export]]
double find_alpha(const int p, const double target) {
  const int number_of_candidates = 300;
  const double lower = 0.0001;
  const double upper = 30;
  double step;
  step = (upper - lower) / static_cast<double>(number_of_candidates - 1);
  double best_x = lower;
  double smallest_difference = std::numeric_limits<double>::infinity();
  for (int i = 0; i < number_of_candidates; ++i) {
    double candidate_x = lower + static_cast<double>(i) * step;
    double loss = 0.0;
    for (int j = 1; j <= p; ++j) {
      loss += candidate_x /
        (static_cast<double>(j) - 1.0 + candidate_x);
    }
    double difference = std::abs(loss - target);
    if (difference < smallest_difference) {
      smallest_difference = difference;
      best_x = candidate_x;
    }
  }
  return best_x;
}
