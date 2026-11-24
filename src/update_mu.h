#ifndef UPDATE_MU_H
#define UPDATE_MU_H

#include <RcppEigen.h>

double update_mu(const Eigen::MatrixXd & X,
                 const Eigen::VectorXd & y,
                 const Eigen::VectorXd & beta,
                 const double sigma2);

#endif  // UPDATE_MU_H
