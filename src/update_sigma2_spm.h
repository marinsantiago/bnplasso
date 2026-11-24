#ifndef UPDATE_SIGMA2_SPM_H
#define UPDATE_SIGMA2_SPM_H

#include <RcppEigen.h>

double update_sigma2_spm_ind_prior(const Eigen::VectorXd & y,
                                   const Eigen::VectorXd & beta);

double update_sigma2_spm_conj_prior(const Eigen::VectorXd & y,
                                    const Eigen::VectorXd & beta,
                                    const Eigen::VectorXd & tau2);

#endif  // UPDATE_SIGMA2_LM_H
