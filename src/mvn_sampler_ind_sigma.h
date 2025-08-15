#ifndef MVN_SAMPLER_IND_SIGMA_H
#define MVN_SAMPLER_IND_SIGMA_H

#include <RcppEigen.h>

Eigen::VectorXd sample_beta_ind_sigma(const Eigen::MatrixXd & X,
                                      const Eigen::MatrixXd & tX,
                                      const Eigen::MatrixXd & tXX,
                                      const Eigen::VectorXd & y,
                                      const Eigen::VectorXd & tau2,
                                      const double sigma2,
                                      const double mu);

#endif  // MVN_SAMPLER_IND_SIGMA_H
