#ifndef UPDATE_BETA_LM_H
#define UPDATE_BETA_LM_H

#include <RcppEigen.h>

Eigen::VectorXd update_beta_lm_ind_sigma(const Eigen::MatrixXd & Xt,
                                         const Eigen::MatrixXd & XtX,
                                         const Eigen::VectorXd & y,
                                         const Eigen::VectorXd & tau2,
                                         const double sigma2,
                                         const double mu);

Eigen::VectorXd update_beta_lm_conj_sigma(const Eigen::MatrixXd & Xt,
                                          const Eigen::MatrixXd & XtX,
                                          const Eigen::VectorXd & y,
                                          const Eigen::VectorXd & tau2,
                                          const double sigma2,
                                          const double mu);

#endif  // UPDATE_BETA_LM_H
