# bnplasso 0.2.2
  - The package now includes a dedicated routine, `bnplasso.spm()`, which 
    implements the sparse normal-means problem. The routine returns an object of
    S3 class, `"spmBayes"`, supported by various class-specific methods. 
  - Other internal functionality improvements.
  - Latest version available at GitHub.

# bnplasso 0.2.1
  - The package now includes the functions `psis.loo()` and `widely.aic()`, 
    which compute the Pareto-smoothed importance sampling leave-one-out 
    information criterion (PSIS-LOO) and the Watanabe–Akaike information 
    criterion (WAIC), respectively.
  - Archived GitHub version: [marinsantiago/bnplasso@6adfa6d]()
  
# bnplasso 0.2.0
  - The functions `blasso.lm()` and `balasso.lm()` have been merged into the
    function, `bnplasso.lm()`. A new argument, `prior`, has been 
    introduced in `bnplasso.lm()`, which specifies the type of shrinkage prior 
    that should be employed. The options are: (1) a nonparametric Bayesian 
    Lasso prior (default), (2) a Bayesian Lasso prior, or (3) a Bayesian 
    adaptive Lasso prior. 
  - The `bnplasso.lm()` function now supports single-precision floating-point 
    calculations for certain internal routines via the `float` argument. By 
    default, the function still uses double point precision.
  - By default, the `bnplasso.lm()` function now includes an intercept term in 
    the regression function.
  - The `bnplasso.lm()` function now returns the log-likelihood of each 
    observation at each MCMC iteration.
  - If some user-supplied hyperparameters are not provided, the `bnplasso.lm()` 
    function will now attempt to automatically determine appropriate values 
    for those hyperparameters.
  - The package now includes the function `get.partition()`, which recovers the 
    partition of the regression coefficients induced by the nonparametric
    Bayesian Lasso.
  - The package now includes the functions `coclust.probs()` and 
    `coclust.point()`, which compute and visualize the matrices of co-clustering
    probabilities and co-clustering point estimates, respectively. 
  - Other internal functionality improvements, including a better handling of 
    numerical instabilities and memory management. 
  - Archived GitHub version: [marinsantiago/bnplasso@4cf966b]()
    
# bnplasso 0.1.0
  - Initial release with core functionality.
  - Implementation as described in 
    <https://doi.org/10.1080/10618600.2025.2572327>
  - Archived GitHub version: [marinsantiago/bnplasso@3c87169]()
