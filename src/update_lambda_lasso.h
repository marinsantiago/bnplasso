#ifndef UPDATE_LAMBDA_LASSO_H
#define UPDATE_LAMBDA_LASSO_H

#include <RcppEigen.h>

Eigen::VectorXd update_lambda2_Bayesian_lasso(
    const Eigen::VectorXd & tau2,
    const double a, 
    const double b
);

Eigen::VectorXd update_lambda2_Bayesian_adapt_lasso(
    const Eigen::VectorXd & tau2,
    const double a,
    const double b
);

#endif  // UPDATE_LAMBDA_LASSO_H
