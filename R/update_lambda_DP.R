# ------------------------------------------------------------------------------
# Update lambda with a Dirichlet process mixture
# ------------------------------------------------------------------------------

updateDP <- function(tau2, lambda2_k, clust.indicators, a, b, alpha) {
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
      lambda2_k[indicator.j] <- rgamma(1, shape = a, rate = b) + 1e-8
    }
    # Update cluster indicators
    clust.indicators[j] <- indicator.j
  }
  # Update cluster membership indicators ---------------------------------------
  #for (j in seq_len(p)) {
  #  K.current <- length(lambda2_k) # Total number of current clusters
  #  # Number of taus in each cluster (excluding tau_j)
  #  n.cluster_minus.j <- tabulate(clust.indicators[-j], nbins = K.current)
  #  probs.j <- c(
  #    # Probabilities that tau_j joins one of the existing clusters
  #    dexp(tau2[j], rate = lambda2_k/2) * n.cluster_minus.j,
  #    # Probability that tau_j joins a new cluster
  #    alpha * (a * b^(a)) / (2 * (b + tau2[j]/2)^(a + 1))
  #  )
  #  probs.j <- probs.j / sum(probs.j) # Normalize the probabilities
  #  # Sample the j-th cluster indicator
  #  indicator.j <- int_sampling(K.current + 1, 1, probs.j)
  #  # If tau_j joins a new cluster, update lambda2 by sampling it from G_0
  #  if (indicator.j == (K.current + 1)) {
  #    lambda2_k[indicator.j] <- rgamma(1, shape = a, rate = b) + 1e-8
  #  }
  #  # Update cluster indicators
  #  clust.indicators[j] <- indicator.j
  #}
  
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
    lambda2_k[k] <- rgamma(1, shape = post_shape, rate = post_rate) + 1e-8
  }
  # If the current cluster is empty, sample lambda2_k from G_0
  lambda2_k[idx_empty_clusters] <- rgamma(
    length(idx_empty_clusters), shape = a, rate = b
  ) + 1e-8
  #K.current <- length(lambda2_k)
  #seq_p <- seq_len(p)
  #for (k in seq_len(K.current)) {
  #  # Indices of tau's that belong to cluster k
  #  idx <- seq_p[clust.indicators == k]
  #  n_k <- length(idx) # Number of taus in the cluster k
  #  if (!n_k) {
  #    # If the current cluster is empty, sample lambda2 from G_0
  #    lambda2_k[k] <- rgamma(1, shape = a, rate = b) + 1e-8
  #  } else {
  #    # Sample from the posterior
  #    post_shape <- a + n_k
  #    post_rate <- b + sum(tau2[idx])/2
  #    lambda2_k[k] <- rgamma(1, shape = post_shape, rate = post_rate) + 1e-8
  #  }
  #}
  unique.clusts <- unique(clust.indicators)
  lambda2_k[lambda2_k <= 1e-8] <- 1e-8
  list(
    clust.indicators = clust.indicators,
    lambda2_k = lambda2_k,
    K.clust = length(unique.clusts),
    unique.clusts = unique.clusts
  )
}


#updateDP <- function(tau2, lambda2_k, clust.indicators, a, b, alpha) {
#  p <- length(tau2)
#  probs.new.clust <- alpha * (a * b^(a)) / (2 * (b + tau2/2)^(a + 1))
#  K.current <- length(lambda2_k) # Current number of current clusters
#  n.cluster <- tabulate(clust.indicators, nbins = K.current) # Clusters size
#  # Update cluster membership indicators --------------------------------------
#  for (j in seq_len(p)) {
#    # Clusters size (excluding tau_j)
#    n.clust_j <- n.cluster
#    n.clust_j[clust.indicators[j]] <- n.clust_j[clust.indicators[j]] - 1
#    # Find the current non-empty clusters
#    non_empty_clusts_minus.j <- seq_len(K.current)[n.clust_j != 0]
#    # Size of the non-empty clusters
#    n.non_empty_clusts_minus.j <- n.clust_j[non_empty_clusts_minus.j]
#    # Lambdas of the non-empty clusters
#    lambda2_k.non_empty <- lambda2_k[non_empty_clusts_minus.j]
#    # Probabilities of joining a non-empty cluster or a new cluster
#    # Recall that the probability of joining an empty cluster is 0
#    probs.j <- c(
#      # Probabilities that tau_j joins one of the existing non-empty clusters
#      dexp(tau2[j], rate = lambda2_k.non_empty/2) * n.non_empty_clusts_minus.j,
#      # Probability that tau_j joins a new cluster
#      probs.new.clust[j]
#    )
#    probs.j <- probs.j / sum(probs.j) # Normalize the probabilities
#    # Sample the j-th cluster indicator (K.current + 1 denotes a new cluster)
#    potential_clusts <- c(non_empty_clusts_minus.j, K.current + 1)
#    K.potnetial <- length(potential_clusts)
#    indicator.j <- potential_clusts[int_sampling(K.potnetial, 1, probs.j)]
#    # Update K.current
#    K.current <- max(K.current, indicator.j)
#    # If tau_j joins a new cluster:
#    if (indicator.j == (K.current + 1)) {
#      # Sample lambda2 from G_0
#      lambda2_k[indicator.j] <- rgamma(1, shape = a, rate = b) + 1e-8
#      # Update n.cluster: tau_j joins the cluster K + 1
#      n.cluster[indicator.j] <-  1
#    } else {
#      # Update n.cluster: tau_j "jumps" from clust.indicators[j] to indicator.j
#      n.cluster[clust.indicators[j]] <- n.cluster[clust.indicators[j]] - 1
#      n.cluster[indicator.j] <- n.cluster[indicator.j] + 1
#    }
#    # Update cluster indicators
#    clust.indicators[j] <- indicator.j
#  }
#  
#  # Update lambda2_k vector ----------------------------------------------------
#  K.current <- max(clust.indicators)
#  non_empty_clusters <- tabulate(clust.indicators, K.current) != 0
#  idx_non_empty_clusters <- which(non_empty_clusters)
#  idx_empty_clusters <- which(!non_empty_clusters)
#  # Loop only over non-empty clusters
#  for (k in idx_non_empty_clusters) {
#    # Sample the corresponding lambda2_k from the full conditional distribution
#    idx <- clust.indicators == k # Which tau's belong to the current cluster
#    n_k <- sum(idx)
#    post_shape <- a + n_k
#    post_rate <- b + sum(tau2[idx])/2
#    lambda2_k[k] <- rgamma(1, shape = post_shape, rate = post_rate) + 1e-8
#  }
#  # If the current cluster is empty, ignore lambda2_k
#  lambda2_k[idx_empty_clusters] <- rep(1e-8, length(idx_empty_clusters))
#  unique.clusts <- unique(clust.indicators)
#  lambda2_k[lambda2_k <= 1e-8] <- 1e-8
#  list(
#    clust.indicators = clust.indicators, lambda2_k = lambda2_k,
#    K.clust = length(unique.clusts), unique.clusts = unique.clusts
#  )
#}