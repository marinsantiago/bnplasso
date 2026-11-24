#ifndef UTILS_H
#define UTILS_H

#include <RcppEigen.h>

bool within_integer_vec(const int x, const std::vector<int> v);

Eigen::MatrixXd center(const Eigen::MatrixXd & x);

Eigen::MatrixXd make_posdef(const Eigen::MatrixXd & m);

double find_alpha(const int p, const double target);

#endif  // UTILS_H
