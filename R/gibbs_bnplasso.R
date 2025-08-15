# ------------------------------------------------------------------------------
# Gibbs sampler the for the nonparametric Bayesian Lasso
# ------------------------------------------------------------------------------

gibbs_bnpLasso_ind_sigma <- function(X, y, a, b, alpha, intercept,
                                     variance.prior.type,
                                     max.iters, burn.in, thin) {
  
  # Get dimensions of the data -------------------------------------------------
  dims <- dim(X)
  n <- dims[1]
  p <- dims[2]
  
  # Initialize data structures to store posterior draws ------------------------
  dd <- max.iters - burn.in # Number of posterior draws to store after burn in
  K.out <- sigma2.out <- rep(NA, dd)
  if (intercept) mu.out <- rep(NA, dd)
  Clust.idx.out <- matrix(NA, nrow = dd, ncol = p)
  lambda2.out <- beta.out <- tau2.out <- Clust.idx.out
  
  # Pre-compute constants ------------------------------------------------------
  tXX <- make_posdef(get_tXX(X)) # cross-product t(X) %*% X
  tX <- t(X)
  
  # Initialize model parameters ------------------------------------------------
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
  tau2 <- abs(rexp(p, rate = mix.rates / 2)) + 1e-06
  sigma2 <- 1
  mu <- 0
  lambda2 <- rep(NA, p)
  
  # MCMC algorithm -------------------------------------------------------------
  pb <- progress::progress_bar$new(
    format = " Gibbs sampling: bnplasso [:bar] :percent in :elapsed",
    total = max.iters, clear = FALSE,
  )
  start <- Sys.time() # Start wall-clock time
  for (iter in seq_len(max.iters)) {
    pb$tick()
    # Update lambda2 (Dirichlet process)
    out.DP <- updateDP(tau2, lambda2_k, clust.indicators, a, b, alpha)
    clust.indicators <- out.DP$clust.indicators
    K.clust <- out.DP$K.clust
    unique.clusts <- out.DP$unique.clusts
    lambda2_k <- out.DP$lambda2_k
    for (clust.k in unique.clusts) {
      lambda2[clust.indicators == clust.k] <- lambda2_k[clust.k]
    }
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
      K.out[current.save] <- K.clust
      sigma2.out[current.save] <- sigma2 
      Clust.idx.out[current.save, ] <- clust.indicators
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
  Clust.idx.out <- Clust.idx.out[thin.seq, ]
  K.out <- K.out[thin.seq]
  if (intercept) mu.out <- mu.out[thin.seq]
  
  # Returns --------------------------------------------------------------------
  out <- list(
    Post.beta = beta.out, Post.sigma2 = sigma2.out, Post.tau2 = tau2.out,
    Post.lambda2 = lambda2.out, Post.clust_idx = Clust.idx.out, Post.K = K.out
  )
  if (intercept) out[["Post.mu"]] <- mu.out
  out[["elapsed"]] <- elapsed
  out
}
