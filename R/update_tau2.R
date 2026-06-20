# ------------------------------------------------------------------------------
# Update the latent parameters, tau2
# ------------------------------------------------------------------------------

# Sample from the full conditional distribution of tau2, assuming an 
# independent variance prior.
sample_tau2_ind_sigma <- function(beta, lambda2) {
  p <- length(beta)
  xi <- sqrt(pmax(abs(lambda2 / (beta^(2))), 1e-08))
  gamm <- actuar::rinvgauss(p, mean = xi, shape = lambda2)
  pmax(abs(1 / gamm), 1e-08)
}

# Sample from the full conditional distribution of sigma2, assuming a 
# conjugate variance prior.
sample_tau2_conj_sigma <- function(beta, lambda2, sigma2) {
  p <- length(beta)
  xi <- sqrt(pmax(abs((lambda2 * sigma2) / (beta^(2))), 1e-08))
  gamm <- actuar::rinvgauss(p, mean = xi, shape = lambda2)
  pmax(abs(1 / gamm), 1e-08)
}
