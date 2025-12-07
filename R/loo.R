# ------------------------------------------------------------------------------
# Leave-one-out and WAIC 
# ------------------------------------------------------------------------------

#' PSIS-LOO
#'
#' This function computes the Pareto-smoothed importance sampling leave-one-out 
#' information criterion (PSIS-LOO) for an object of class \code{lmBayes}. 
#' The PSIS-LOO is converted into deviance scale so that it is comparable 
#' with other information criteria like AIC and BIC.
#' 
#' @param object An object of class \code{'lmBayes'}.
#' @param comp.r_eff Logical. Whether the function should estimate the relative 
#'    effective sample size for the likelihood of each observation. If 
#'    \code{FALSE}, all the relative effective sample sizes are set to one.
#'    Default is \code{TRUE}.
#' 
#' @return The estimated PSIS-LOO.
#' 
#' @author Santiago Marin
#'
psis.loo <- function(object, comp.r_eff = TRUE) {
  if (!is.lmBayes(object)) stop("object should be of class 'lmBayes'")
  if (not.logic(comp.r_eff)) stop("comp.r_eff must be logical of length one")
  if (!comp.r_eff) {
    return(suppressWarnings(loo::loo(object$loglik))$estimates[3,1])
  }
  likelihood <- exp(object$loglik)
  rel_eff <- loo::relative_eff(likelihood, chain_id = rep(1, nrow(likelihood)))
  suppressWarnings(loo::loo(object$loglik, r_eff = rel_eff))$estimates[3,1]
}


#' WAIC
#'
#' This function computes the Watanabe–Akaike information criterion (WAIC)
#' for an object of class \code{lmBayes}. The WAIC is converted into deviance 
#' scale so that it is comparable with other information criteria 
#' like AIC and BIC.
#' 
#' @param object An object of class \code{'lmBayes'}.
#' 
#' @return The estimated WAIC.
#' 
#' @author Santiago Marin
#'
widely.aic <- function(object) {
  if (!is.lmBayes(object)) stop("object should be of class 'lmBayes'")
  suppressWarnings(loo::waic(object$loglik))$estimates[3,1]
}
