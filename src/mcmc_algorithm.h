#ifndef MCMC_ALGORITHM_H
#define MCMC_ALGORITHM_H

#include <RcppEigen.h>

Rcpp::List mcmc_sampler(const Eigen::MatrixXd & X,
                        const Eigen::VectorXd & y,
                        const double a,
                        const double b,
                        const double alpha,
                        const bool intercept,
                        const std::string shrinkage_prior, 
                        const std::string variance_prior_type,
                        const bool polya_urn_sampler,
                        const bool sp_means,
                        const int max_iters,
                        const int burn_in,
                        const int thin,
                        const bool verbose);

#endif  // MCMC_ALGORITHM_H
