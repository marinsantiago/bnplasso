# ------------------------------------------------------------------------------
# Gibbs sampler that implements either:
#   (1) Nonparametric Bayesian Lasso.
#   (2) Bayesian Lasso.
#   (3) Bayesian Adaptive Lasso.
# ------------------------------------------------------------------------------

mcmc_sampler <- function(X, y, a, b, alpha, intercept, penalty.type, 
                         variance.prior.type, max.iters, burn.in, thin, 
                         polya, float, sp.means) {
  
  # Get dimensions of the data -------------------------------------------------
  if (!sp.means) {
    dims <- dim(X)
    n <- dims[1]
    p <- dims[2]
  } else {
    n <- p <- length(y)
    intercept <- FALSE
  }
  
  # Initialize data structures to store posterior draws ------------------------
  # Number of draws after burn-in and thinning
  dd <- floor((max.iters - burn.in) / thin)
  K.out <- sigma2.out <- rep(NA, dd)
  if (intercept) mu.out <- rep(NA, dd)
  Clust.idx.out <- matrix(NA, nrow = dd, ncol = p)
  lambda2.out <- beta.out <- tau2.out <- Clust.idx.out
  if (!polya) {
    # Blocked Gibbs sampler
    K <- 1e+02
    nu.out <- matrix(NA, nrow = dd, ncol = K - 1)
    omega.out <- matrix(NA, nrow = dd, ncol = K)
  }
  
  # Pre-compute constants ------------------------------------------------------
  if (!sp.means) {
    tXX <- make_posdef(get_tXX(X)) # cross-product t(X) %*% X
    tX <- t(X)
  }
  
  # Initialize model parameters ------------------------------------------------
  if (penalty.type == "bnp.lasso") {
    if (polya) {
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
      rm(mix.rates)
    } else { # Blocked Gibbs sampler
      # Sample all lambda2 candidates from G_0
      lambda2.all <- rgamma(K, shape = a, rate = b)
      # Generate nu and the stick-breaking weights (omega)
      nu <- rbeta(K - 1, 1, alpha)
      nu.full <- c(nu, 1) # The last nu is always one
      omega <- rep(NA, K)
      omega[1] <- nu.full[1]
      for (k in 2:K) omega[k] <- nu.full[k] * prod(1 - nu.full[1:(k-1)])
      # Sample the cluster allocations given the above omega
      clust.indicators <- sample.int(K, p, replace = T, prob = omega)
      # Sample tau2
      tau2 <- pmax(rexp(p, rate = lambda2.all[clust.indicators] / 2), 1e-06)
      rm(nu.full)
    }
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
      if (polya) {
        # Generalized Polya urn sampling scheme
        out.DP.Polya <- updateDP.Polya(
          tau2, lambda2_k, clust.indicators, a, b, alpha
        )
        clust.indicators <- out.DP.Polya$clust.indicators
        K.clust <- out.DP.Polya$K.clust
        unique.clusts <- out.DP.Polya$unique.clusts
        lambda2_k <- out.DP.Polya$lambda2_k
        for (clust.k in unique.clusts) {
          lambda2[clust.indicators == clust.k] <- lambda2_k[clust.k]
        }
      } else {
        # Blocked Gibbs sampling scheme
        out.DP.Gibbs <- updateDP.blockedGibbs(
          tau2 = tau2, lambda2.all = lambda2.all, 
          clust.indicators = clust.indicators, nu = nu, omega = omega,
          a = a, b = b, alpha = alpha, K = K, float = float
        )
        clust.indicators <- out.DP.Gibbs$clust.indicators
        K.clust <- out.DP.Gibbs$K.clust
        unique.clusts <- out.DP.Gibbs$unique.clusts
        lambda2.all <- out.DP.Gibbs$lambda2.all
        nu <- out.DP.Gibbs$nu
        omega <- out.DP.Gibbs$omega
        lambda2 <- out.DP.Gibbs$lambda2
      }
    } else if (penalty.type == "b.lasso") {
      lambda2 <- sample_lambda2_b.lasso(a, b, tau2)
    } else if (penalty.type == "b.adapt.lasso") {
      lambda2 <- sample_lambda2_ba.lasso(a, b, tau2)
    }
    # Update beta, tau2, and sigma2, assuming an independent variance prior
    if (variance.prior.type == "independent") {
      # Linear regression:
      if (!sp.means) {
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
      } else { # Sparse means problem:
        # Update betas
        betas <- sample_beta_ind_sigma_sp_means(y, tau2, sigma2) 
        # Update tau2
        tau2 <- sample_tau2_ind_sigma(betas, lambda2)
        # Update sigma2
        sigma2 <- sample_sigma2_ind_prior_sp_means(y, betas) 
      }
    } else if (variance.prior.type == "conjugate") {
      # Linear regression:
      if (!sp.means) {
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
      } else { # Sparse means problem:
        # Update betas
        betas <- sample_beta_conj_sigma_sp_means(y, tau2, sigma2)
        # Update tau2
        tau2 <- sample_tau2_conj_sigma(betas, lambda2, sigma2)
        # Update sigma2
        sigma2 <- sample_sigma2_conj_prior_sp_means(y, betas, tau2)
      }
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
          if (!polya) {
            nu.out[current.save, ] <- nu
            omega.out[current.save, ] <- omega
          }
        }
      }
    }
    if (iter %% 1000 == 0) gc() # Collect cache
  }
  end <- Sys.time()
  elapsed <- difftime(end, start, units = "secs")
  
  # Returns --------------------------------------------------------------------
  out <- list(
    post.beta = beta.out, post.sigma2 = sigma2.out, 
    post.tau2 = tau2.out, post.lambda2 = lambda2.out
  )
  if (intercept) out[["post.mu"]] <- mu.out
  if (penalty.type == "bnp.lasso") {
    out[["post.K"]] <- K.out
    out[["post.clust_idx"]] <- Clust.idx.out
    if (!polya) {
      out[["post.nu"]] <- nu.out
      out[["post.omega"]] <- omega.out
    }
  }
  out[["elapsed"]] <- elapsed
  out
}
