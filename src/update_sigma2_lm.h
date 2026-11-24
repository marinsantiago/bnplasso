#ifndef UPDATE_SIGMA2_LM_H
#define UPDATE_SIGMA2_LM_H

#include <RcppEigen.h>

double update_sigma2_ind_prior(const Eigen::MatrixXd & X,
                               const Eigen::VectorXd & y,
                               const Eigen::VectorXd & beta,
                               const double mu,
                               const bool intercept);

double update_sigma2_conj_prior(const Eigen::MatrixXd & X,
                                const Eigen::VectorXd & y,
                                const Eigen::VectorXd & beta,
                                const Eigen::VectorXd & tau2,
                                const double mu,
                                const bool intercept);

#endif  // UPDATE_SIGMA2_LM_H
