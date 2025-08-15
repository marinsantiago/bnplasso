#' Expected log pointwise predictive density (elppd)
#'
#' This function computes the expected log pointwise predictive density (elppd) 
#' for an object of class \code{lmBayes}. Given a test (held-out) set, 
#' \eqn{\{y_{\text{test},i},\,
#' \mathbf{x}_{\text{test},i}\}_{i=1}^{n_{\text{test}}}}, the elppd is defined 
#' as 
#' \deqn{\text{elppd} = \frac{1}{n_{\text{test}}}\sum_{i=1}^{n_{\text{test}}}
#' \log\left(\frac{1}{S}\sum_{s=1}^{S}p\left({y}_{\text{test}, i} | 
#' \mathbf{x}_{\text{test}, i}^{_{'}}\boldsymbol{\beta}_{(s)}, 
#' \sigma^{2}_{(s)}\right)\right),}
#' where \eqn{p()} is the likelihood of the observed data and 
#' \eqn{\boldsymbol{\beta}_{(s)}} and \eqn{\sigma^{2}_{(s)}} are the 
#' \eqn{s}-th draws of the coefficient vector and the sampling variance in the 
#' MCMC algorithm, respectively. 
#'
#' @param object An object of class \code{lmBayes}.
#' @param X.new A matrix of predictors from the test data, where each row is an 
#'   observation. 
#' @param y.new A vector of responses from the test data of size 
#'   \eqn{n_{\text{test}}}.
#' 
#' @return The elppd as a numeric scalar.
#'
#' @references
#'
#' A. Gelman, J. Carlin, H. Stern, D. Dunson, and A. Vehtari (2013), 
#' \emph{Bayesian Data Analysis}, 3. Chapman & Hall/CRC.
#' 
#' @author Santiago Marin
#' 
elppd <- function(object, X.new, y.new) {
  
  # Input validation -----------------------------------------------------------
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  n.draws <- object$n.draws
  n.preds <- object$n.preds
  chck <- (!is.numeric(X.new)) || (!is.matrix(X.new))
  if(chck) stop("X.new must be a numeric matrix.")
  X.new.dims <- dim(X.new)
  if(X.new.dims[2] != n.preds) stop("X.new has an incorrect number of columns.")
  chck <- (!is.numeric(y.new)) || (!is.vector(y.new))
  if (chck) stop("y.new must be a numeric vector.")
  n.test <- length(y.new)
  chck <- X.new.dims[1] != n.test
  mssg <- "Number of observations in y.new does not match the number of rows"
  if(chck) stop(paste(mssg, "in X.new."))
  
  # Compute elppd --------------------------------------------------------------
  
  # Pre-compute all linear predictors (without the intercept) for the new data.
  linPreds <- tcrossprod(X.new, object$Post.beta)
  # Add the intercept if needed
  if (object$intercept) {
    linPreds <- linPreds + matrix(
      data = object$Post.mu, nrow = nrow(linPreds),
      ncol = ncol(linPreds), byrow = TRUE
    )
  }
  post_sigmas <- sqrt(abs(object$Post.sigma2) + 1e-16)
  # log pointwise predictive density
  lppd <- sapply(
    seq_len(n.test), \(i) {
      log(mean(dnorm(y.new[i], linPreds[i,], post_sigmas) + 1e-25))
    }
  )
  elppd.out <- mean(lppd)
  rm(lppd, linPreds, post_sigmas); gc()
  elppd.out
}
