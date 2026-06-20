# ------------------------------------------------------------------------------
# Routines for co-clustering analysis
# ------------------------------------------------------------------------------

#' Posterior co-clustering probabilities plot
#'
#' This function plots the posterior co-clustering probabilities of the 
#' regression coefficients
#' 
#' @param S.mcmc A matrix of size \code{n.draws}-by-\code{n.preds}, where each 
#'    row indicates to which cluster the regression coefficients belong to.
#' @param return.mat logical. Whether the matrix of posterior co-clustering 
#'    probabilities should be returned. Default is \code{FALSE}
#' @param viridis.pal Character describing which "viridis" palette should be
#'    employed. Options are "A", "B", "C", "D", "E", "F", "G", and "H".
#' @param main Title of the plot.
#' @param axis.idx Sequence of number to appear in the axis of the plot.
#' @param xy.labs A title for the x and y axes.
#' 
#' @author Santiago Marin
#' 
coclust.probs <- function(S.mcmc, return.mat = FALSE, viridis.pal = "C", 
                          main = NULL, axis.idx = NULL, xy.labs = NULL) {
  
  # Input validation -----------------------------------------------------------
  if (not.mat(S.mcmc))  stop("S.mcmc must be a valid numeric matrix")
  if (not.logic(return.mat)) stop("return.mat must be logical")
  if (!(viridis.pal %in% LETTERS[1:7])) {
    stop(paste("viridis.pal must be", paste(LETTERS[1:7], collapse = ", ")))
  }
  
  # Pre-compute constants ------------------------------------------------------
  S.dims <- dim(S.mcmc) 
  n.draws <- S.dims[1]
  n.preds <- S.dims[2]
  
  # Compute co-clustering matrix -----------------------------------------------
  coclustering <- matrix(0, nrow = n.preds, ncol = n.preds)
  for (iter in seq_len(n.draws)) {
    cluster <- S.mcmc[iter, ]
    coclustering <- coclustering + outer(cluster, cluster, FUN = "==")
  }
  coclustering.out <- coclustering / n.draws
  diag(coclustering.out) <- rep(1, n.preds) # The diagonal is always one
  if (return.mat) return(coclustering.out)
  
  # Generate plot --------------------------------------------------------------
  idxs <- if (is.null(axis.idx)) {
    floor(seq(1, n.preds, length.out = min(5L, n.preds)))
  } else { axis.idx } 
  main.title <- if (is.null(main)) {
    "Posterior co-clustering probabilities"
  } else { main }
  xy.labs <- if (is.null(xy.labs)) "Coefficient index" else xy.labs
  # Compute: "df <- reshape2::melt(coclustering.out)" in pure base R!
  expand.idxs <- arrayInd(
    seq_along(coclustering.out), .dim = dim(coclustering.out)
  )
  df <- data.frame(
    Var1 = expand.idxs[,1],
    Var2 = expand.idxs[,2],
    probability = as.vector(coclustering.out)
  )
  ggplot(df, aes(Var1, Var2, fill = probability)) +
    geom_tile() +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    ggtitle(main.title) + 
    labs(x = xy.labs, y = xy.labs) +
    scale_x_continuous(breaks = idxs, expand = c(0, 0)) +
    scale_y_continuous(breaks = idxs, expand = c(0, 0)) +
    viridis::scale_fill_viridis(name = "Prob.", option = viridis.pal) +
    theme(axis.title.x = element_text(margin = margin(t = 8))) +
    theme(axis.title.y = element_text(margin = margin(r = 8))) +
    theme(axis.text = element_text(size = 10)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    theme(legend.title = element_text(margin = margin(b = 15)))
}


#' Posterior co-clustering plot
#'
#' Plots the posterior co-clustering of the regression coefficients
#' 
#' @param S A vector of size \code{n.preds}, indicating to which cluster the 
#'    corresponding regression coefficients belong to.
#' @param return.mat logical. Whether the matrix of posterior co-clustering 
#'    probabilities should be returned. Default is \code{FALSE}
#' @param viridis.pal Character describing which "viridis" palette should be
#'    employed. Options are "A", "B", "C", "D", "E", "F", "G", and "H".
#' @param main Title of the plot.
#' @param axis.idx Sequence of number to appear in the axis of the plot.
#' @param xy.labs A title for the x and y axes.
#' 
#' @author Santiago Marin
#' 
coclust.point <- function(S, return.mat = FALSE, viridis.pal = "C", 
                          main = NULL, axis.idx = NULL, xy.labs = NULL) {
  
  # Input validation -----------------------------------------------------------
  if (not.logic(return.mat)) stop("return.mat must be logical")
  if (not.y(S))  stop("S must be a valid numeric vector")
  if (!(viridis.pal %in% LETTERS[1:7])) {
    stop(paste("viridis.pal must be", paste(LETTERS[1:7], collapse = ", ")))
  }
  
  # Pre-compute constants ------------------------------------------------------
  n.preds <- length(S)
  
  # Compute co-clustering matrix -----------------------------------------------
  coclusts <- outer(S, S, FUN = "==") * 1
  diag(coclusts) <- rep(1, n.preds)  # The diagonal is always one
  if (return.mat) return(coclusts)
  
  # Gen. plot ------------------------------------------------------------------
  lower.col <- viridis::viridis(9, option = viridis.pal)[1]
  upper.col <- viridis::viridis(9, option = viridis.pal)[9]
  main.title <- if (is.null(main)) "Co-clustering" else main
  xy.labs <- if (is.null(xy.labs)) "Coefficient index" else xy.labs
  if (is.null(axis.idx)) {
    idx <- rep("", n.preds)
    idx.seq <- floor(seq(1, n.preds, length.out = 5))
    idx[idx.seq] <- as.character(idx.seq)
  } else {
    idx <- rep("", n.preds)
    idx[axis.idx] <- as.character(axis.idx)
  }
  # Compute: "df <- reshape2::melt(coclusts)" in pure base R!
  expand.idxs <- arrayInd(seq_along(coclusts), .dim = dim(coclusts))
  df <- data.frame(
    Var1 = expand.idxs[,1], Var2 = expand.idxs[,2], value = as.vector(coclusts)
  )
  df$`co-cluster` <- factor(
    ifelse(df$value == 1, "yes", "no"), 
    levels = c("yes", "no")
  )
  ggplot(df, aes(Var1, Var2, fill = `co-cluster`)) + 
    geom_tile() +
    scale_fill_manual(values = c("yes" = upper.col, "no" = lower.col)) + 
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    ggtitle(main.title) + 
    labs(x = xy.labs, y = xy.labs) +
    scale_x_discrete(limits = factor(1:n.preds), labels = idx, breaks = idx) + 
    scale_y_discrete(limits = factor(1:n.preds), labels = idx, breaks = idx) + 
    theme(axis.title.x = element_text(margin = margin(t = 8))) +
    theme(axis.title.y = element_text(margin = margin(r = 8))) +
    theme(axis.text = element_text(size = 10)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    theme(legend.title = element_text(margin = margin(b = 15)))
}
