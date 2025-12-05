# ------------------------------------------------------------------------------
# Gibbs sampler that implements either:
#   (1) Nonparametric Bayesian Lasso.
#   (2) Bayesian Lasso.
#   (3) Bayesian Adaptive Lasso.
# ------------------------------------------------------------------------------

mcmc_sampler <- function(X, y, a, b, alpha, intercept, penalty.type,
                         variance.prior.type, max.iters, burn.in, thin, float) {
  
  # Get dimensions of the data -------------------------------------------------
  dims <- dim(X)
  n <- dims[1]
  p <- dims[2]
  
  # Initialize data structures to store posterior draws ------------------------
  # Number of draws after burn-in and thinning
  dd <- floor((max.iters - burn.in) / thin) 
  K.out <- sigma2.out <- rep(NA, dd)
  if (intercept) mu.out <- rep(NA, dd)
  Clust.idx.out <- matrix(NA, nrow = dd, ncol = p)
  lambda2.out <- beta.out <- tau2.out <- Clust.idx.out
  
  # Pre-compute constants ------------------------------------------------------
  tXX <- make_posdef(get_tXX(X)) # cross-product t(X) %*% X
  tX <- t(X)
  
  # Initialize model parameters ------------------------------------------------
  if (penalty.type == "bnp.lasso") {
    # Randomly assign coefficients to some initial clusters
    n.clust.init <- 4L
    clust.indicators <- sample.int(n.clust.init, p, replace = T)
    # Consider a sequence of lambda2_k values applying different shrinkage
    lambda2_k <- c(0.01, 0.1, 0.5, 200)
    # Get the corresponding exponential mixing rates
    mix.rates <- rep(NA, p)
    for (cl in seq_len(n.clust.init)) {
      mix.rates[clust.indicators == cl] <- lambda2_k[cl]
    }
    # Draw each tau_j from their respective prior distributions
    tau2 <- pmax(abs(rexp(p, rate = mix.rates / 2)), 1e-06)
    lambda2 <- rep(NA, p)
  } else if (penalty.type == "b.lasso") {
    # Sample lambda2 from its respective prior distribution
    lambda2 <- rep(max(abs(rgamma(1, shape = a, rate = b)), 1e-08), p)
    # Draw each tau_j from their respective prior distributions
    tau2 <- pmax(abs(rexp(p, rate = lambda2 / 2)), 1e-06)
  } else if (penalty.type == "b.adapt.lasso") {
    # Sample each lambda2 from their respective prior distributions
    lambda2 <- pmax(abs(rgamma(p, shape = a, rate = b)), 1e-08)
    # Draw each tau_j from their respective prior distributions
    tau2 <- pmax(abs(rexp(p, rate = lambda2 / 2)), 1e-06)
  }
  sigma2 <- 1
  mu <- 0
  
  # MCMC algorithm -------------------------------------------------------------
  pb <- progress::progress_bar$new(
    format = " Gibbs sampler: [:bar] :percent in :elapsed",
    total = max.iters, clear = FALSE,
  )
  start <- Sys.time() # Start wall-clock time
  for (iter in seq_len(max.iters)) {
    pb$tick()
    # Update lambda2
    if (penalty.type == "bnp.lasso") {
      out.DP <- updateDP(tau2, lambda2_k, clust.indicators, a, b, alpha)
      clust.indicators <- out.DP$clust.indicators
      K.clust <- out.DP$K.clust
      unique.clusts <- out.DP$unique.clusts
      lambda2_k <- out.DP$lambda2_k
      for (clust.k in unique.clusts) {
        lambda2[clust.indicators == clust.k] <- lambda2_k[clust.k]
      }
    } else if (penalty.type == "b.lasso") {
      lambda2 <- sample_lambda2_b.lasso(a, b, tau2)
    } else if (penalty.type == "b.adapt.lasso") {
      lambda2 <- sample_lambda2_ba.lasso(a, b, tau2)
    }
    # Update beta, tau2, and sigma2, assuming an independent variance prior
    if (variance.prior.type == "independent") {
      # Update betas
      if (float) {
        betas <- sample_beta_ind_sigma_float(X, tX, tXX, y, tau2, sigma2, mu)
      } else {
        betas <- sample_beta_ind_sigma(X, tX, tXX, y, tau2, sigma2, mu)
      }
      # Update tau2
      tau2 <- sample_tau2_ind_sigma(betas, lambda2)
      # Update sigma2
      sigma2 <- sample_sigma2_ind_prior(X, y, betas, mu, intercept)
    } else if (variance.prior.type == "conjugate") {
      # Update betas
      if (float) {
        betas <- sample_beta_conj_sigma_float(X, tX, tXX, y, tau2, sigma2, mu)
      } else {
        betas <- sample_beta_conj_sigma(X, tX, tXX, y, tau2, sigma2, mu)
      }
      # Update tau2
      tau2 <- sample_tau2_conj_sigma(betas, lambda2, sigma2)
      # Update sigma2
      sigma2 <- sample_sigma2_conj_prior(X, y, betas, tau2, mu, intercept)
    }
    # Update mu
    mu <- if (intercept) sample_mu(X, y, betas, sigma2) else 0
    # Store draws after burn-in and thinning
    if (iter > burn.in) {
      if (iter %% thin == 0) {
        current.save <- iter - burn.in
        sigma2.out[current.save] <- sigma2 
        lambda2.out[current.save, ] <- lambda2
        beta.out[current.save, ] <- betas
        tau2.out[current.save, ] <- tau2
        if (intercept) mu.out[current.save] <- mu
        if (penalty.type == "bnp.lasso") {
          K.out[current.save] <- K.clust
          Clust.idx.out[current.save, ] <- clust.indicators
        }
      }
    }
    if (iter %% 1000 == 0) gc() # Collect cache
  }
  end <- Sys.time()
  elapsed <- end - start
  
  # Returns --------------------------------------------------------------------
  out <- list(
    Post.beta = beta.out, Post.sigma2 = sigma2.out, 
    Post.tau2 = tau2.out, Post.lambda2 = lambda2.out
  )
  if (intercept) out[["Post.mu"]] <- mu.out
  if (penalty.type == "bnp.lasso") {
    out[["Post.K"]] <- K.out
    out[["Post.clust_idx"]] <- Clust.idx.out
  }
  out[["elapsed"]] <- elapsed
  out
}
