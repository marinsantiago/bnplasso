# ------------------------------------------------------------------------------
# Update the latent parameters, tau2
# ------------------------------------------------------------------------------

# Sample from the full conditional distribution of tau2, assuming an 
# independent variance prior.
sample_tau2_ind_sigma <- function(beta, lambda2) {
  p <- length(beta)
  xi <- sqrt(abs(lambda2 / (beta^(2))) + 1e-8)
  gamm <- actuar::rinvgauss(p, mean = xi, shape = lambda2)
  tau2 <- abs(1 / gamm) + 1e-8
  #tau2[tau2 <= 1e-8] <- 1e-8
  tau2
}


# Sample from the full conditional distribution of sigma2, assuming a 
# conjugate variance prior.
sample_tau2_conj_sigma <- function(beta, lambda2, sigma2) {
  p <- length(beta)
  xi <- sqrt(abs((lambda2 * sigma2) / (beta^(2))) + 1e-8)
  gamm <- actuar::rinvgauss(p, mean = xi, shape = lambda2)
  tau2 <- abs(1 / gamm) + 1e-8
  #tau2[tau2 <= 1e-8] <- 1e-8
  tau2
}
