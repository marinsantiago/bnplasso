# ------------------------------------------------------------------------------
# Update the intercept term, mu
# ------------------------------------------------------------------------------

# Sample from the full conditional distribution of mu
sample_mu <- function(X, y, beta, sigma2) {
  n <- length(y)
  post.mean <- mean(y - X %*% beta)
  post.sd <- sqrt(abs(sigma2 / n) + 1e-8)
  rnorm(1, post.mean, post.sd)
}
