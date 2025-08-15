# ------------------------------------------------------------------------------
# Gibbs sampler the for the Bayesian Lasso
# ------------------------------------------------------------------------------

gibbs_bLasso <- function(X, y, a, b, intercept,
                         variance.prior.type,
                         max.iters, burn.in, thin) {
  
  # Get dimensions of the data -------------------------------------------------
  dims <- dim(X)
  n <- dims[1]
  p <- dims[2]
  
  # Initialize data structures to store posterior draws ------------------------
  dd <- max.iters - burn.in # Number of posterior draws to store after burn in
  sigma2.out <- rep(NA, dd)
  if (intercept) mu.out <- rep(NA, dd)
  lambda2.out <- beta.out <- tau2.out <- matrix(NA, nrow = dd, ncol = p)
  
  # Pre-compute constants ------------------------------------------------------
  tXX <- make_posdef(get_tXX(X)) # cross-product t(X) %*% X
  tX <- t(X)
  
  # Initialize model parameters ------------------------------------------------
  # Sample lambda2 from its respective prior distribution
  lambda2 <- rep(abs(rgamma(1, shape = a, rate = b)) + 1e-8, p)
  # Draw each tau_j from their respective prior distributions
  tau2 <- abs(rexp(p, rate = lambda2 / 2)) + 1e-06
  sigma2 <- 1
  mu <- 0
  
  # MCMC algorithm -------------------------------------------------------------
  pb <- progress::progress_bar$new(
    format = " Gibbs sampling: blasso [:bar] :percent in :elapsed",
    total = max.iters, clear = FALSE,
  )
  start <- Sys.time() # Start wall-clock time
  for (iter in seq_len(max.iters)) {
    pb$tick()
    # Update lambda2 (Bayesian Lasso)
    lambda2 <- sample_lambda2_b.lasso(a, b, tau2)
    # Assuming an independent variance prior
    if (variance.prior.type == "independent") {
      # Update betas
      betas <- sample_beta_ind_sigma(X, tX, tXX, y, tau2, sigma2, mu)
      # Update tau2
      tau2 <- sample_tau2_ind_sigma(betas, lambda2)
      # Update sigma2
      sigma2 <- sample_sigma2_ind_prior(X, y, betas, mu, intercept)
    } else if (variance.prior.type == "conjugate") {
      # Update betas
      betas <- sample_beta_conj_sigma(X, tX, tXX, y, tau2, sigma2, mu)
      # Update tau2
      tau2 <- sample_tau2_conj_sigma(betas, lambda2, sigma2)
      # Update sigma2
      sigma2 <- sample_sigma2_conj_prior(X, y, betas, tau2, mu, intercept)
    }
    # Update mu
    mu <- if (intercept) sample_mu(X, y, betas, sigma2) else 0
    
    # Store posterior draws
    if (iter > burn.in) {
      current.save <- iter - burn.in
      sigma2.out[current.save] <- sigma2 
      lambda2.out[current.save, ] <- lambda2
      beta.out[current.save, ] <- betas
      tau2.out[current.save, ] <- tau2
      if (intercept) mu.out[current.save] <- mu
    }
    if (iter %% 1000 == 0) gc()
  }
  end <- Sys.time()
  elapsed <- end - start
  
  # Apply thinning -------------------------------------------------------------
  thin.seq <- seq(from = 1, to = dd, by = thin)
  beta.out <- beta.out[thin.seq, ]
  sigma2.out <- sigma2.out[thin.seq]
  tau2.out <- tau2.out[thin.seq, ]
  lambda2.out <- lambda2.out[thin.seq, ]
  if (intercept) mu.out <- mu.out[thin.seq]
  
  # Returns --------------------------------------------------------------------
  out <- list(
    Post.beta = beta.out, Post.sigma2 = sigma2.out,
    Post.tau2 = tau2.out, Post.lambda2 = lambda2.out
  )
  if (intercept) out[["Post.mu"]] <- mu.out
  out[["elapsed"]] <- elapsed
  out
}
