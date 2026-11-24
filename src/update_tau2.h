#ifndef UPDATE_TAU2_LM_H
#define UPDATE_TAU2_LM_H

#include <RcppEigen.h>

Eigen::VectorXd update_tau2_ind_sigma(const Eigen::VectorXd & beta,
                                      const Eigen::VectorXd & lambda2);

Eigen::VectorXd update_tau2_conj_sigma(const Eigen::VectorXd & beta,
                                       const Eigen::VectorXd & lambda2,
                                       const double sigma2);

#endif  // UPDATE_TAU2_LM_H
