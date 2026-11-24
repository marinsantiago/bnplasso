#ifndef UPDATE_BETA_SPM_H
#define UPDATE_BETA_SPM_H

#include <RcppEigen.h>

Eigen::VectorXd update_beta_spm_ind_sigma(const Eigen::VectorXd & y,
                                          const Eigen::VectorXd & tau2,
                                          const double sigma2);

Eigen::VectorXd update_beta_spm_conj_sigma(const Eigen::VectorXd & y,
                                           const Eigen::VectorXd & tau2,
                                           const double sigma2);

#endif  // UPDATE_BETA_SPM_H
