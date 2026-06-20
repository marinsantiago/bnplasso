# ------------------------------------------------------------------------------
# Parameter updates for the sparse means problems
# ------------------------------------------------------------------------------

# Independent sigma case -------------------------------------------------------
sample_beta_ind_sigma_sp_means <- function(y, tau2, sigma2) {
  p <- length(tau2)
  # Compute Phi^(-1) (only store the diagonal as a vector)
  tau2.inv <- pmax(1 / tau2, 1e-15)
  sigma2.inv <- max(1 / sigma2, 1e-15)
  Phi.inv <- pmax(1 / (sigma2.inv + tau2.inv), 1e-15)
  # Compute the mean
  m <- Phi.inv * y / sigma2
  # Generate random draws
  rnorm(p, mean = m, sd = sqrt(Phi.inv))
}

sample_sigma2_ind_prior_sp_means <- function(y, beta) {
  n <- length(y)
  ssr <- crossprod(y - beta)
  scale.s2 <- ssr / 2
  shape.s2 <- (n) / 2
  out <- actuar::rinvgamma(1, shape = shape.s2, scale = scale.s2)
  max(out, 1e-08)
}

# Conjugate sigma case ---------------------------------------------------------
sample_beta_conj_sigma_sp_means <- function(y, tau2, sigma2) {
  p <- length(tau2)
  # Compute Phi^(-1) (only store the diagonal as a vector)
  tau2.inv <- pmax(1 / tau2, 1e-15)
  Phi.inv <- pmax(1 / (1 + (tau2.inv)), 1e-15)
  # Compute the standard deviation
  sd <- sqrt(pmax(sigma2 * Phi.inv, 1e-15))
  # Compute the mean
  m <- Phi.inv * y
  # Generate random draws
  rnorm(p, mean = m, sd = sd)
}

sample_sigma2_conj_prior_sp_means <- function(y, beta, tau2) {
  n <- length(y)
  p <- length(beta)
  # Quadratic form of betas
  tau2.inv <- pmax(1 / tau2, 1e-08)
  quad_form_beta <- crossprod(beta * tau2.inv, beta)
  ssr <- crossprod(y - beta)
  scale.s2 <- (ssr + quad_form_beta) / 2
  shape.s2 <- (n + p) / 2
  out <- actuar::rinvgamma(1, shape = shape.s2, scale = scale.s2)
  max(out, 1e-08)
}
