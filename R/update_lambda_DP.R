# ------------------------------------------------------------------------------
# Update lambda with a Dirichlet process mixture
# ------------------------------------------------------------------------------

# Update the DP assuming a Polya urn scheme
updateDP.Polya <- function(tau2, lambda2_k, clust.indicators, a, b, alpha) {
  p <- length(tau2)
  # Update cluster membership indicators ---------------------------------------
  for (j in seq_len(p)) {
    K.current <- length(lambda2_k) # Total number of current clusters
    # Number of taus in each cluster (excluding tau_j)
    n.cluster_minus.j <- tabulate(clust.indicators[-j], nbins = K.current)
    # Compute the probabilities of the non-empty clusters only or a new cluster
    non_empty_clusts_minus.j <- which(n.cluster_minus.j != 0) 
    n.non_empty_clusts_minus.j <- n.cluster_minus.j[non_empty_clusts_minus.j]
    lambda2_k.non_empty <- lambda2_k[non_empty_clusts_minus.j]
    probs.j <- c(
      # Probabilities that tau_j joins one of the existing non-empty clusters
      dexp(tau2[j], rate = lambda2_k.non_empty/2) * n.non_empty_clusts_minus.j,
      # Probability that tau_j joins a new cluster
      alpha * (a * b^(a)) / (2 * (b + tau2[j]/2)^(a + 1))
    )
    probs.j <- probs.j / sum(probs.j) # Normalize the probabilities
    # Sample the j-th cluster indicator (K.current + 1 denotes a new cluster)
    potential_clusts <- c(non_empty_clusts_minus.j, K.current + 1)
    K.potnetial <- length(potential_clusts)
    indicator.j <- potential_clusts[int_sampling(K.potnetial, 1, probs.j)]
    # If tau_j joins a new cluster, update lambda2 by sampling it from G_0
    if (indicator.j == (K.current + 1)) {
      lambda2_k[indicator.j] <- max(rgamma(1, shape = a, rate = b), 1e-08)
    }
    # Update cluster indicators
    clust.indicators[j] <- indicator.j
  }

  # Update lambda2_k vector ----------------------------------------------------
  K.current <- max(clust.indicators)
  non_empty_clusters <- tabulate(clust.indicators, K.current) != 0
  idx_non_empty_clusters <- which(non_empty_clusters)
  idx_empty_clusters <- which(!non_empty_clusters)
  # Loop only over non-empty clusters
  for (k in idx_non_empty_clusters) {
    # Sample the corresponding lambda2_k from the full conditional distribution
    idx <- clust.indicators == k # Which tau's belong to the current cluster
    n_k <- sum(idx)
    post_shape <- a + n_k
    post_rate <- b + sum(tau2[idx])/2
    lambda2_k[k] <- max(rgamma(1, shape = post_shape, rate = post_rate), 1e-08)
  }
  # If the current cluster is empty, sample lambda2_k from G_0
  lambda2_k[idx_empty_clusters] <- pmax(
    rgamma(length(idx_empty_clusters), shape = a, rate = b), 1e-08
  )
  unique.clusts <- unique(clust.indicators)
  #lambda2_k[lambda2_k <= 1e-08] <- 1e-08
  lambda2_k <- pmax(lambda2_k, 1e-08)
  list(
    clust.indicators = clust.indicators, lambda2_k = lambda2_k,
    K.clust = length(unique.clusts), unique.clusts = unique.clusts
  )
}

# Update the DP assuming a blocked Gibbs sampler scheme
updateDP.blockedGibbs <- function(tau2, lambda2.all, clust.indicators, 
                                  nu, omega, a, b, alpha, K, float) {
  p <- length(tau2)
  # Step 1: Update all the lambdas ---------------------------------------------
  allocated.comps <- sort(unique(clust.indicators))
  nonallocated.comps <- setdiff(1:K, allocated.comps)
  K.non <- length(nonallocated.comps)
  K.yes <- length(allocated.comps)
  # Sample the lambdas associated with non-allocated components from G_0
  lambda2.all[nonallocated.comps] <- rgamma(K.non, shape = a, rate = b)
  # Updated lambdas from allocated components
  for (k in allocated.comps) {
    idx.k <- clust.indicators == k
    p.k <- sum(idx.k)
    sum.tau2.k <- sum(tau2[idx.k])
    pos.shape <- a + p.k
    pos.rate <- b + 0.5 * sum.tau2.k
    lambda2.all[k] <- max(rgamma(1, shape = pos.shape, rate = pos.rate), 1e-08)
  }
  
  # Step 2: Update cluster allocations -----------------------------------------
  #for (j in seq_len(p)) {
  #  prob.j <- omega * dexp(tau2[j], rate = lambda2.all / 2) 
  #  prob.j <- prob.j / sum(prob.j)
  #  clust.indicators[j] <- int_sampling(K, 1, prob.j)
  #}
  clust.indicators <- if (!float) {
    step2_blockedGibbs_dbl(tau2, lambda2.all, omega, K)
  } else {
    step2_blockedGibbs_flt(tau2, lambda2.all, omega, K)
  } 

  # Step 3: Update stick-breaking weights --------------------------------------
  nu <- rep(NA, K - 1)
  for (k in seq_len(K - 1)) {
    idx.k <- clust.indicators == k
    idx.kplus <- clust.indicators %in% (k + 1):K
    p.k <- sum(idx.k)
    p.kplus <- sum(idx.kplus)
    nu[k] <- rbeta(1, 1 + p.k, alpha + p.kplus)
  }
  nu.full <- c(nu, 1)
  omega <- rep(NA, K)
  omega[1] <- nu.full[1]
  for (k in 2:K) omega[k] <- nu.full[k] * prod(1 - nu.full[1:(k-1)])
  # Prepare the returns
  lambda2 <- rep(NA, p)
  unique.clusts <- unique(clust.indicators)
  for (k in unique.clusts) lambda2[clust.indicators == k] <- lambda2.all[k]
  list(
    nu = nu, omega = omega, lambda2 = lambda2, 
    K.clust = length(unique.clusts), unique.clusts = unique.clusts, 
    clust.indicators = clust.indicators, lambda2.all = lambda2.all
  )
}
