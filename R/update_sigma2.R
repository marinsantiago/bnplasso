# ------------------------------------------------------------------------------
# Update the sampling variance, sigma2
# ------------------------------------------------------------------------------

# Sample from the full conditional distribution of sigma2, assuming an 
# independent variance prior. "intercept" is a logical variable 
# indicating if an intercept should be included in the model.
sample_sigma2_ind_prior <- function(X, y, beta, mu, intercept) {
  n <- length(y)
  if (intercept) {
    res <- y - mu - X %*% beta
    ssr <- crossprod(res)
    scale.s2 <- ssr / 2
    shape.s2 <- n / 2
  } else {
    res <- y - X %*% beta
    ssr <- crossprod(res)
    scale.s2 <- ssr / 2
    shape.s2 <- (n - 1) / 2
  }
  out <- actuar::rinvgamma(1, shape = shape.s2, scale = scale.s2)
  if (out <= 1e-08) return(1e-08) 
  out
}

# Sample from the full conditional distribution of sigma2, assuming a 
# conjugate variance prior. "intercept" is a logical variable 
# indicating if an intercept should be included in the model.
sample_sigma2_conj_prior <- function(X, y, beta, tau2, mu, intercept) {
  n <- length(y)
  p <- length(beta)
  # Quadratic form of betas
  itau2 <- abs(1 / tau2) + 1e-8
  beta_itau2 <- beta * itau2
  quad_form_beta <- crossprod(beta_itau2, beta)
  if (intercept) {
    res <- y - mu - X %*% beta
    ssr <- crossprod(res)
    scale.s2 <- (ssr + quad_form_beta) / 2
    shape.s2 <- (n + p) / 2
  } else {
    res <- y - X %*% beta
    ssr <- crossprod(res)
    scale.s2 <- (ssr + quad_form_beta) / 2
    shape.s2 <- (n + p - 1) / 2
  }
  out <- actuar::rinvgamma(1, shape = shape.s2, scale = scale.s2)
  if (out <= 1e-08) return(1e-08) 
  out
}
