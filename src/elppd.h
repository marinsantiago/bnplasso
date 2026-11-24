#ifndef ELPPD_H
#define ELPPD_H

#include <RcppEigen.h>

double elppd(const Rcpp::List & object, 
             const Eigen::MatrixXd & X_new,
             const Eigen::VectorXd & y_new);

#endif  // ELPPD_H
