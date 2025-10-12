#' Fit linear regression models with a Bayesian Lasso prior
#'
#' This function fits linear regression models with a Bayesian Lasso prior as in
#' Park and Casella (2008).
#'
#' @param X A matrix of predictors of dimension \eqn{n}-by-\eqn{p}, where each 
#'   of the \eqn{n} rows is an observation vector.
#' @param y Response variable. It should be a numeric vector size \eqn{n}.
#' @param a A positive scalar corresponding to the \bold{shape} parameter in the 
#'   gamma distribution used as hyper-prior in the shrinkage parameter,
#'   \eqn{\lambda^{2}}.
#' @param b A positive scalar corresponding to the \bold{rate} parameter in the 
#'   gamma distribution used as hyper-prior in the shrinkage parameter,
#'   \eqn{\lambda^{2}}.
#' @param intercept Logical. If \code{TRUE}, an intercept term is included in the
#'   model; otherwise, the intercept is integrated out. If \code{TRUE}, the 
#'   prior on the intercept would be a non-informative prior of the form 
#'   \eqn{p(\mu)\propto 1}. Default is \code{FALSE}.
#' @param variance.prior.type A character string denoting whether the variance 
#'   on the sampling variance should be an independent-type prior or a 
#'   conjugate-type prior. See Moran et al. (2019) for details. The options are
#'    either \code{"independent"} or \code{"conjugate"}. 
#'    Default is \code{"conjugate"} as in Park and Casella (2008).
#' @param max.iters A positive integer corresponding to the total number of 
#'   MCMC iterations. Default is 6000.  
#' @param burn.in A positive integer corresponding to the number of draws 
#'   discarded as burn-in. It should be smaller than \code{max.iters}.
#'   Default is 1000.
#' @param thin A positive integer specifying the period for saving samples.
#'   Default is 1.
#'
#' @return An object of S3 class, "lmBayes", containing:
#' \itemize{
#'   \item \code{Post.beta}: A matrix of size \code{n.draws}-by-\code{n.preds}, 
#'   where each row is a posterior draw of the regression coefficients.
#'   \item \code{Post.sigma2}: A vector of size \code{n.draws}, where each 
#'   element is a posterior draw of the sampling variance. 
#'   \item \code{Post.tau2}: A matrix of size \code{n.draws}-by-\code{n.preds}, 
#'   where each row is a posterior draw of the latent parameters 
#'   \eqn{\tau_{j}^{2}}.
#'   \item \code{Post.lambda2}: A matrix of size 
#'   \code{n.draws}-by-\code{n.preds}, where each row is a posterior draw of the
#'   shrinkage parameters \eqn{\lambda_{j}^{2}}.
#'   \item \code{Post.mu}: A vector of size \code{n.draws}, where each element 
#'   is a posterior draw of the intercept term.
#'   \item \code{elapsed}: The elapsed (wall-clock) time of the MCMC sampler.
#'   \item \code{a}: The \bold{shape} parameter in the gamma distribution 
#'   used as a hyper-prior.
#'   \item \code{b}: The \bold{rate} parameter in the gamma distribution 
#'   used as a hyper-prior.  
#'   \item \code{intercept}: Whether or not an intercept term was included in 
#'   the model.
#'   \item \code{variance.prior.type}: Whether the variance on the sampling 
#'   variance was an independent-type prior or a conjugate-type prior.
#'   \item \code{max.iters}: The total number of MCMC iterations.
#'   \item \code{burn.in}: The number of draws discarded as burn-in.
#'   \item \code{thin}: The period for saving draws.
#'   \item \code{n.obs}: The sample size.   
#'   \item \code{n.preds}: The number of predictors.
#'   \item \code{n.draws}: The number of posterior draws after burn-in and
#'   thinning.
#'   \item \code{X}: Matrix of predictors.
#'   \item \code{y}: Vector of responses.
#'   \item \code{post.pred.fitted.values}: A matrix of size 
#'   \code{n.draws}-by-\code{n.obs}, where each row is a draw from the posterior 
#'   predictive distribution of the fitted values.
#'   \item \code{post.pred.residuals}: A matrix of size 
#'   \code{n.draws}-by-\code{n.obs}, where each row is a draw from the posterior 
#'   predictive distribution of the residuals.  
#' }
#' 
#' @references
#' 
#' G. E. Moran, V. Rockova, and E. I. George (2019), Variance Prior Forms for 
#' High-Dimensional Bayesian Variable Selection. \emph{Bayesian Analysis},
#' 14(4):1091-1119.
#' 
#' T. Park, and G. Casella (2008), The Bayesian Lasso. 
#' \emph{Journal of the American Statistical Association}, 103(482):681-686.
#' 
#' @author Santiago Marin
#' 
blasso.lm <- function(X, y, a, b,
                      intercept = FALSE,
                      variance.prior.type = "conjugate",
                      max.iters = 6000L, burn.in = 1000L, thin = 1L) {
  
  # Input validation -----------------------------------------------------------
  if (!is.numeric(X) || !is.matrix(X)) stop("X must be a numeric matrix.")
  if (!is.numeric(y) || !is.vector(y)) stop("y must be a numeric vector.")
  chck <- !(length(y) == nrow(X))
  if (chck) stop("Length of y does not match the number of rows in X.")
  chck.not.scalar <- \(x)  (!is.numeric(x)) || (length(x) != 1) || (x <= 0)
  if (chck.not.scalar(a)) stop("a must be a positive scalar.")
  if (chck.not.scalar(b)) stop("b must be a positive scalar.")
  chck <- (!is.logical(intercept)) || (length(intercept) != 1)
  if (chck) stop("intercept must be logical of size 1.")
  chck <- !(variance.prior.type %in% c("independent", "conjugate"))
  if (chck) stop("variance.prior.type must be 'independent' or 'conjugate'.")
  chck.not.int <- \(x) (x %% 1) != 0 || (x <= 0) || (length(x) != 1)
  if (chck.not.int(max.iters)) stop("max.iters must be a positive integer.")
  chck <- chck.not.int(burn.in) || max.iters <= burn.in
  if (chck) stop("burn.in must be a positive integer smaller than max.iters.")
  if (chck.not.int(thin)) stop("thin must be a positive integer.")
  rm(chck, chck.not.scalar, chck.not.int); gc()

  # Run Gibbs sampler ----------------------------------------------------------
  blasso_out <- gibbs_bLasso(
    X = X, y = y, a = a, b = b, intercept = intercept,
    variance.prior.type = variance.prior.type,
    max.iters = max.iters, burn.in = burn.in, thin = thin
  )
  
  # Prepare the returns --------------------------------------------------------
  dimsX <- dim(X)
  n.obs <- dimsX[1]
  n.preds <- dimsX[2] 
  n.draws <- nrow(blasso_out$Post.beta)
  blasso_out[["a"]] <- a
  blasso_out[["b"]] <- b
  blasso_out[["intercept"]] <- intercept
  blasso_out[["variance.prior.type"]] <- variance.prior.type
  blasso_out[["max.iters"]] <- max.iters
  blasso_out[["burn.in"]] <- burn.in
  blasso_out[["thin"]] <- thin
  blasso_out[["n.obs"]] <- n.obs
  blasso_out[["n.preds"]] <- n.preds
  blasso_out[["n.draws"]] <- n.draws # After burn-in and thinning 
  
  # Posterior predictive fitted values and residuals ---------------------------
  # Pre-compute all linear predictors (without the intercept)
  linPreds <- tcrossprod(X, blasso_out$Post.beta)
  # Add the intercept if needed
  if (intercept) {
    linPreds <- linPreds + matrix(
      data = blasso_out$Post.mu, nrow = nrow(linPreds),
      ncol = ncol(linPreds), byrow = TRUE
    )
  }
  # Note: In "post_pred_fits" and "post_pred_res", each row corresponds to an
  # MCMC draw and each column to an observation.
  post_pred_fits <- matrix(data = NA, nrow = n.draws, ncol = n.obs)
  post_pred_res <- post_pred_fits
  post_sigmas <- sqrt(abs(blasso_out$Post.sigma2) + 1e-16)
  for (s in seq_len(n.draws)) {
    fit_val_s <- linPreds[,s] + rnorm(n.obs, 0, post_sigmas[s])
    post_pred_fits[s,] <- fit_val_s
    post_pred_res[s,] <- y - fit_val_s
  }
  rm(linPreds, post_sigmas, fit_val_s); gc()
  blasso_out[["X"]] <- X
  blasso_out[["y"]] <- y
  blasso_out[["post.pred.fitted.values"]] <- post_pred_fits
  blasso_out[["post.pred.residuals"]] <- post_pred_res
  class(blasso_out) <- "lmBayes"
  blasso_out
}
