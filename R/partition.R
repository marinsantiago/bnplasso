# ------------------------------------------------------------------------------
# Get the partition of the regression coefficients induced by the "bnplasso"
# ------------------------------------------------------------------------------

#' Recover a partition
#'
#' This function recovers the partition of the regression coefficients induced 
#' by the nonparametric Bayesian Lasso. This function is based on the 
#' implementation from the \code{"BNPmix"} package. 
#' 
#' @param S.mcmc A matrix of size \code{n.draws}-by-\code{n.preds}, where each 
#'    row indicates to which cluster the regression coefficients belong to.
#' @param loss A loss function defined on the space of partitions. It can be 
#'    either the variation of information loss function (\code{"VI"}) or the
#'    Binder loss function (\code{"Binder"}). Default is \code{"VI"}. See 
#'    Wade and Ghahramani (2018) for additional details.
#' 
#' @references
#' 
#' S. Wade and Z. Ghahramani (2018), Bayesian cluster analysis: Point 
#' estimation and credible balls (with discussion). \emph{Bayesian Analysis},
#' 13(2):559-626.
#' 
#' @author Santiago Marin
#' 
get.partition <- function(S.mcmc, loss = "VI") {
  
  # Input validation -----------------------------------------------------------
  if (not.mat(S.mcmc))  stop("S.mcmc must be a valid numeric matrix")
  
  cls.draws <- apply(S.mcmc, 1, \(x) as.numeric(as.factor(x))) |> t()
  #object <- list(clust = cls.draws - 1)
  object <- list(clust = cls.draws)
  class(object) <- "BNPdens"
  BNPmix.out <- BNPmix::partition(object, dist = loss)  
  min.lower.bound.expected.loss <- which.min(BNPmix.out$scores)
  BNPmix.out$partitions[min.lower.bound.expected.loss,]
}
