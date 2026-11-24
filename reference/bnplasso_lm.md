# Fit linear regression models with a nonparametric Bayesian Lasso prior

This function fits linear regression models with a nonparametric
Bayesian Lasso prior as in Marin et al. (2026).

## Usage

``` r
bnplasso_lm(
  X,
  y,
  intercept = TRUE,
  shrinkage_prior = NULL,
  a = NULL,
  b = NULL,
  alpha = NULL,
  variance_prior_type = NULL,
  max_iters = 10000L,
  burn_in = 2000L,
  thin = 2L,
  polya_urn_sampler = TRUE,
  verbose = TRUE
)
```

## Arguments

- X:

  A numeric matrix of predictors of dimension \\n\\-by-\\p\\, where each
  of the \\n\\ rows is an observation vector.

- y:

  Response variable. A numeric vector size \\n\\.

- intercept:

  Logical. If `TRUE`, an intercept term is included in the model;
  otherwise, the intercept is integrated out. If `TRUE`, the prior on
  the intercept would be a non-informative prior of the form
  \\p(\mu)\propto 1\\. Default is `TRUE`.

- shrinkage_prior:

  A character string denoting which type of shrinkage prior should be
  employed. The options are either: (1) `"bnplasso"` which implements
  the nonparametric Bayesian Lasso as in Marin et al. (2026), (2)
  `"b_lasso"` which implements the Bayesian Lasso as in Park and Casella
  (2008), or (3) `"b_adapt_lasso"` which implements the Bayesian
  adaptive Lasso as in Leng et al. (2014). Default is `"bnplasso"`.

- a:

  A positive scalar. If `prior = "bnplasso"`, it corresponds to the
  **shape** parameter in the gamma distribution used as a centering
  measure in the DP prior. If `prior = "b_lasso"`, it corresponds to the
  **shape** parameter in the gamma distribution used as prior for the
  shrinkage parameter, \\\lambda^{2}\\. If `prior = "b_adapt_lasso"`, it
  corresponds to the **shape** parameter in the gamma distribution used
  as prior for the shrinkage parameters, \\\lambda\_{j}^{2}\\, for
  \\j\in\\1,\dots,p\\\\. If not provided, the function will attempt to
  determine an appropriate value.

- b:

  A positive scalar. If `prior = "bnplasso"`, it corresponds to the
  **rate** parameter in the gamma distribution used as a centering
  measure in the DP prior. If `prior = "b_lasso"`, it corresponds to the
  **rate** parameter in the gamma distribution used as prior for the
  shrinkage parameter, \\\lambda^{2}\\. If `prior = "b_adapt_lasso"`, it
  corresponds to the **rate** parameter in the gamma distribution used
  as prior for the shrinkage parameters, \\\lambda\_{j}^{2}\\, for
  \\j\in\\1,\dots,p\\\\. If not provided, the function will attempt to
  determine an appropriate value.

- alpha:

  A positive scalar. If `prior = "bnplasso"`, it corresponds to the
  **concentration** parameter in the DP prior. If `prior = "b_lasso"` or
  `prior = "b_adapt_lasso"`, the argument is ignored. If not provided,
  the function will attempt to determine an appropriate value.

- variance_prior_type:

  A character string denoting whether the prior on the sampling variance
  should be an independent-type prior or a conjugate-type prior. See
  Moran et al. (2019) for details. The options are either
  `"independent"` or `"conjugate"`. If not provided, the function will
  attempt to determine an appropriate prior.

- max_iters:

  A positive integer corresponding to the total number of MCMC
  iterations. Default is 10000.

- burn_in:

  A positive integer corresponding to the number of draws discarded as
  burn-in. It should be smaller than `max_iters`. Default is 20000.

- thin:

  A positive integer specifying the period for saving samples. Default
  is 2.

- polya_urn_sampler:

  logical. Whether a generalized Pólya urn sampling scheme or a blocked
  Gibbs sampling scheme should be employed to update the Dirichlet
  process mixture in the MCMC algorithm. If `TRUE`, a generalized Pólya
  urn sampling scheme would be used. Default is `TRUE`. This argument is
  ignored if `prior != "bnp.lasso"`

- verbose:

  Logical. Whether a progress bar should be printed in the the console.
  Default is `TRUE`.

## Value

An object of S3 class, `'lmBayes'`, containing:

- `post.beta`: A matrix of size `n.draws`-by-`n.preds`, where each row
  is a posterior draw of the regression coefficients.

- `post.sigma2`: A vector of size `n.draws`, where each element is a
  posterior draw of the sampling variance.

- `post.tau2`: A matrix of size `n.draws`-by-`n.preds`, where each row
  is a posterior draw of the latent parameters \\\tau\_{j}^{2}\\.

- `post.lambda2`: A matrix of size `n.draws`-by-`n.preds`, where each
  row is a posterior draw of the shrinkage parameters
  \\\lambda\_{j}^{2}\\.

- `post.clust_idx`: If `prior = "bnp.lasso"`, a matrix of size
  `n.draws`-by-`n.preds`, where each row indicates to which cluster the
  regression coefficients belong to.

- `post.K`: If `prior = "bnp.lasso"`, a vector of size `n.draws`, where
  each element indicates the number of clusters in the corresponding
  MCMC iteration.

- `post.mu`: If `intercept = "TRUE"`, a vector of size `n.draws`, where
  each entry is a posterior draw of the intercept term.

- `elapsed`: The elapsed (wall-clock) time of the MCMC sampler, in
  seconds.

- `a`: If `prior = "bnp_lasso"`, the **shape** parameter in the gamma
  distribution used as a centering measure in the DP prior. If
  `prior = "b_lasso"`, the **shape** parameter in the gamma distribution
  used as prior in the shrinkage parameter, \\\lambda^{2}\\. If
  `prior = "b_adapt_lasso"`, the **shape** parameter in the gamma
  distribution used as prior in the shrinkage parameters,
  \\\lambda\_{j}^{2}\\, for \\j\in\\1,\dots,p\\\\.

- `b`: If `prior = "bnplasso"`, the **rate** parameter in the gamma
  distribution used as a centering measure in the DP prior. If
  `prior = "b_lasso"`, the **rate** parameter in the gamma distribution
  used as prior in the shrinkage parameter, \\\lambda^{2}\\. If
  `prior = "b_adapt_lasso"`, the **rate** parameter in the gamma
  distribution used as prior in the shrinkage parameters,
  \\\lambda\_{j}^{2}\\, for \\j\in\\1,\dots,p\\\\.

- `alpha`: If `prior = "bnplasso"`, the **concentration** parameter in
  the Dirichlet process prior.

- `intercept`: Whether an intercept was included in the model.

- `variance.prior.type`: Whether the variance on the sampling variance
  was an independent-type prior or a conjugate-type prior.

- `max.iters`: The total number of MCMC iterations.

- `burn.in`: The number of draws discarded as burn-in.

- `thin`: The period for saving draws.

- `n.obs`: The sample size.

- `n.preds`: The number of predictors.

- `n.draws`: The number of posterior draws after burn-in and thinning.

- `X`: Matrix of predictors.

- `y`: Vector of responses.

- `loglik`: Matrix of size `n.draws`-by-`n.obs` with the log-likelihood
  of each observation at each MCMC iteration.

- `post.pred.fitted.values`: A matrix of size `n.draws`-by-`n.obs`,
  where each row is a draw from the posterior predictive distribution of
  the fitted values.

- `post.pred.residuals`: A matrix of size `n.draws`-by-`n.obs`, where
  each row is a draw from the posterior predictive distribution of the
  residuals.

## References

C. Leng, MN. Tran, and D. Nott (2014), Bayesian adaptive Lasso. *Ann
Inst Stat Math*, 66:221-244

S. Marin, B. Long,and A. H. Westveld (2026), Adaptive Shrinkage with a
Nonparametric Bayesian Lasso.*Journal of Computational and Graphical
Statistics*, 35(2):854-864. doi:10.1080/10618600.2025.2572327

G. E. Moran, V. Rockova, and E. I. George (2019), Variance Prior Forms
for High-Dimensional Bayesian Variable Selection. *Bayesian Analysis*,
14(4):1091-1119.

T. Park, and G. Casella (2008), The Bayesian Lasso. *Journal of the
American Statistical Association*, 103(482):681-686.

## Author

Santiago Marin
