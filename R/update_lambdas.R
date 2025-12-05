# ------------------------------------------------------------------------------
# Update the shrinkage parameters from the Bayesian Lasso and the
# Bayesian adaptive Lasso
# ------------------------------------------------------------------------------

# Sample from the full conditional distribution of each lambda_{j}^{2} under the
# Bayesian **adaptive** Lasso
sample_lambda2_ba.lasso <- function(a, b, tau2) {
  p <- length(tau2)
  post.shape <- rep(abs(a + 1) + 1e-8, p)
  post.rate <- abs(b + (tau2 / 2)) + 1e-8
  lambda2 <- rgamma(p, shape = post.shape, rate = post.rate)
  pmax(abs(lambda2), 1e-8)
}


# Sample from the full conditional distribution of lambda^{2} under the 
# Bayesian Lasso
sample_lambda2_b.lasso <- function(a, b, tau2) {
  p <- length(tau2)
  post.shape <- abs(a + p) + 1e-8
  post.rate <- abs(b + sum(tau2)/2) + 1e-8
  rnd_draw <- rgamma(1, shape = post.shape, rate = post.rate)
  rnd_draw <- max(abs(rnd_draw), 1e-08)
  rep(rnd_draw, p)
}
