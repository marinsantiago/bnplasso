// -----------------------------------------------------------------------------
// Random number generation from the inverse Gaussian distribution
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include "rinvgauss.h"


// Random number generation from the inverse-Gaussian(mu, lambda) distribution.
// [[Rcpp::export]]
double rinvgauss(const double mu, const double lambda) {
  const double n = R::rnorm(0.0, 1.0);
  const double muy = mu * n * n;
  const double h = 0.5 * mu / lambda;
  // Larger root
  const double x2 = mu + h * (muy + std::sqrt(4.0 * lambda * muy + muy * muy)); 
  // Smaller root
  const double x1 = mu * mu / x2;
  const double u = R::runif(0.0, 1.0);
  if (u * (mu + x1) <= mu) {
    return x1;
  } else {
    return x2;
  }
}
