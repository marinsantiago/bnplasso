#ifndef MVN_SAMPLER_JEFFREYS_SIGMA_H
#define MVN_SAMPLER_JEFFREYS_SIGMA_H

#include <RcppEigen.h>

Eigen::VectorXd sample_beta_conj_sigma(const Eigen::MatrixXd & X,
                                       const Eigen::MatrixXd & tX,
                                       const Eigen::MatrixXd & tXX,
                                       const Eigen::VectorXd & y,
                                       const Eigen::VectorXd & tau2,
                                       const double sigma2,
                                       const double mu);

Eigen::VectorXf sample_beta_conj_sigma_float(const Eigen::MatrixXf & X,
                                             const Eigen::MatrixXf & tX,
                                             const Eigen::MatrixXf & tXX,
                                             const Eigen::VectorXf & y,
                                             const Eigen::VectorXf & tau2,
                                             const float sigma2,
                                             const float mu);

#endif  // MVN_SAMPLER_JEFFREYS_SIGMA_H
