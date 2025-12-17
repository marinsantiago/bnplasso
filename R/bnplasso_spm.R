#' Fit sparse means models with a nonparametric Bayesian Lasso prior
#'
#' This function fits sparse means models with a nonparametric 
#' Bayesian Lasso prior (Marin et al., 2025+).
#'
#' @param y Response variable. It should be a numeric vector size \eqn{n}.
#' @param prior A character string denoting which type of shrinkage prior
#'   should be employed. The options are either: (1) \code{"bnp.lasso"} which
#'   implements the nonparametric Bayesian Lasso as in Marin et al. (2025+), 
#'   (2) \code{"b.lasso"} which implements the Bayesian Lasso as in Park and 
#'   Casella (2008), or (3) \code{"b.adapt.lasso"} which implements the Bayesian
#'   adaptive Lasso as in Leng et al. (2014). Default is \code{"bnp.lasso"}.
#' @param a A positive scalar. If \code{prior = "bnp.lasso"}, it 
#'   corresponds to the \bold{shape} parameter in the gamma distribution used 
#'   as a centering measure in the DP prior. If \code{prior = "b.lasso"}, 
#'   it corresponds to the \bold{shape} parameter in the gamma distribution used 
#'   as prior in the shrinkage parameter, \eqn{\lambda^{2}}. If 
#'   \code{prior = "b.adapt.lasso"}, it corresponds to the \bold{shape}
#'   parameter in the gamma distribution used as prior in the shrinkage 
#'   parameters, \eqn{\lambda_{j}^{2}}, for \eqn{j\in\{1,\dots,p\}}. If not 
#'   provided, the function will attempt to determine an appropriate value.
#' @param b A positive scalar. If \code{prior = "bnp.lasso"}, it 
#'   corresponds to the \bold{rate} parameter in the gamma distribution used 
#'   as a centering measure in the DP prior. If \code{prior = "b.lasso"}, 
#'   it corresponds to the \bold{rate} parameter in the gamma distribution used 
#'   as prior in the shrinkage parameter, \eqn{\lambda^{2}}. If 
#'   \code{prior = "b.adapt.lasso"}, it corresponds to the \bold{rate}
#'   parameter in the gamma distribution used as prior in the shrinkage 
#'   parameters, \eqn{\lambda_{j}^{2}}, for \eqn{j\in\{1,\dots,p\}}. If not 
#'   provided, the function will attempt to determine an appropriate value.
#' @param alpha A positive scalar. If \code{prior = "bnp.lasso"}, it 
#'   corresponds to the \bold{concentration} parameter in the DP prior. If 
#'   \code{prior = "b.lasso"} or \code{prior = "b.adapt.lasso"},
#'   the argument is ignored. If not provided, the function will attempt to 
#'   determine an appropriate value.
#' @param variance.prior.type A character string denoting whether the prior 
#'   on the sampling variance should be an independent-type prior or a 
#'   conjugate-type prior. See Moran et al. (2019) for details. The options are 
#'   either \code{"independent"} or \code{"conjugate"}. If not provided, the 
#'   function will attempt to determine an appropriate prior.
#' @param max.iters A positive integer corresponding to the total number of 
#'   MCMC iterations. Default is 6000.  
#' @param burn.in A positive integer corresponding to the number of draws 
#'   discarded as burn-in. It should be smaller than \code{max.iters}.
#'   Default is 1000.
#' @param thin A positive integer specifying the period for saving samples.
#'   Default is 1.
#' @param polya logical. Whether a generalized Pólya urn sampling scheme or a 
#'   blocked Gibbs sampling scheme should be employed to update the Dirichlet 
#'   process mixture in the MCMC algorithm. If \code{TRUE}, a generalized 
#'   Pólya urn sampling scheme would be used. Default is \code{FALSE}. 
#'   Important: This argument is ignored if \code{prior != "bnp.lasso"}
#' @param float Logical. If \code{TRUE}, some internal routines will use single 
#'   point precision for improved computational efficiency at the expense of 
#'   numerical accuracy. If \code{FALSE}, double point precision will be used. 
#'   Default is \code{FALSE}.
#'
#' @return An object of S3 class, \code{'spmBayes'}, containing:
#' \itemize{
#'   \item \code{post.beta}: A matrix of size \code{n.draws}-by-\code{n.preds}, 
#'   where each row is a posterior draw of the regression coefficients.
#'   \item \code{post.sigma2}: A vector of size \code{n.draws}, where each 
#'   element is a posterior draw of the sampling variance. 
#'   \item \code{post.tau2}: A matrix of size \code{n.draws}-by-\code{n.preds}, 
#'   where each row is a posterior draw of the latent parameters 
#'   \eqn{\tau_{j}^{2}}.
#'   \item \code{post.lambda2}: A matrix of size 
#'   \code{n.draws}-by-\code{n.preds}, where each row is a posterior draw of the
#'   shrinkage parameters \eqn{\lambda_{j}^{2}}.
#'   \item \code{post.clust_idx}: If \code{prior = "bnp.lasso"}, a matrix of 
#'   size \code{n.draws}-by-\code{n.preds}, where each row indicates to which 
#'   cluster the regression coefficients belong to.
#'   \item \code{post.K}: If \code{prior = "bnp.lasso"}, a vector of size 
#'   \code{n.draws}, where each element indicates the number of clusters in the 
#'   corresponding MCMC iteration.
#'   \item \code{elapsed}: The elapsed (wall-clock) time of the MCMC sampler, 
#'   in seconds.  
#'   \item \code{a}: If \code{prior = "bnp.lasso"}, the \bold{shape} parameter 
#'   in the gamma distribution used as a centering measure in the DP prior. 
#'   If \code{prior = "b.lasso"}, the \bold{shape} parameter in the gamma 
#'   distribution used as prior in the shrinkage parameter, \eqn{\lambda^{2}}. 
#'   If \code{prior = "b.adapt.lasso"}, the \bold{shape} parameter in the gamma 
#'   distribution used as prior in the shrinkage parameters, 
#'   \eqn{\lambda_{j}^{2}}, for \eqn{j\in\{1,\dots,p\}}.
#'   \item \code{b}: If \code{prior = "bnp.lasso"}, the \bold{rate} parameter 
#'   in the gamma distribution used as a centering measure in the DP prior. 
#'   If \code{prior = "b.lasso"}, the \bold{rate} parameter in the gamma 
#'   distribution used as prior in the shrinkage parameter, \eqn{\lambda^{2}}. 
#'   If \code{prior = "b.adapt.lasso"}, the \bold{rate} parameter in the gamma 
#'   distribution used as prior in the shrinkage parameters, 
#'   \eqn{\lambda_{j}^{2}}, for \eqn{j\in\{1,\dots,p\}}.
#'   \item \code{alpha}: If \code{prior = "bnp.lasso"}, the \bold{concentration} 
#'   parameter in the Dirichlet process prior.
#'   \item \code{variance.prior.type}: Whether the variance on the sampling 
#'   variance was an independent-type prior or a conjugate-type prior.
#'   \item \code{max.iters}: The total number of MCMC iterations.
#'   \item \code{burn.in}: The number of draws discarded as burn-in.
#'   \item \code{thin}: The period for saving draws.
#'   \item \code{n.obs}: The sample size.   
#'   \item \code{n.preds}: The number of predictors.
#'   \item \code{n.draws}: The number of posterior draws after burn-in and
#'   thinning.
#'   \item \code{y}: Vector of responses.
#'   \item \code{loglik}: Matrix of size \code{n.draws}-by-\code{n.obs} with the 
#'   log-likelihood of each observation at each MCMC iteration.
#'   \item \code{post.pred}: A matrix of size \code{n.draws}-by-\code{n.obs}, 
#'   where each row is a draw from the posterior predictive distribution.
#' }
#' 
#' @references
#'
#' C. Leng, MN. Tran, and D. Nott (2014), Bayesian adaptive Lasso. 
#' \emph{Ann Inst Stat Math}, 66:221-244
#'
#' S. Marin, B. Long,and A. H. Westveld (2025+), Adaptive Shrinkage with a 
#' Nonparametric Bayesian Lasso.\emph{Journal of Computational and Graphical 
#' Statistics}. doi:10.1080/10618600.2025.2572327
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
#' @encoding UTF-8
#'
bnplasso.spm <- function(y, prior = NULL, a = NULL, b = NULL, 
                         alpha = NULL, variance.prior.type = NULL, 
                         max.iters = 6000L, burn.in = 1000L, thin = 1L,
                         polya = FALSE, float = FALSE) {
  
  # Input validation -----------------------------------------------------------
  if(not.y(y)) stop("y is not a valid response vector")
  if (!is.null(prior)) {if (not.pr(prior)) stop("not a valid shrinakge prior")}
  if (!is.null(variance.prior.type)) {
    if (not.var(variance.prior.type)) stop("not a valid variance.prior.type")
  }
  for (s in c("a", "b", "alpha")) {
    if (!is.null(get(s))) {
      if (not.ps(get(s))) stop(paste(s, " must be a positive scalar"))
    }
  }
  for (i in c("max.iters", "burn.in", "thin")) {
    if (not.int(get(i))) stop(paste(i, " must be a positive integer"))
  }
  if (max.iters < burn.in) stop("burn.in must be smaller than max.iters")
  if (not.logic(float)) stop("float must be logical of length one")
  if (not.logic(polya)) stop("polya must be logical of length one")
  gc()
  
  # Set default hyper-parameters (if not provided) -----------------------------
  n <- p <- length(y)
  if (is.null(prior)) prior <- "bnp.lasso"
  if (is.null(a)) a <- 0.1
  if (is.null(b)) b <- if (prior == "bnp.lasso") max(0.01 / n, 1e-10) else 0.01
  if (is.null(alpha)) {
    if (prior == "bnp.lasso") {
      # Choose alpha so that prior K is approx. 3
      candidate.alphas <- seq(0.0001, 30, length.out = 300)
      prior.K <- sapply(candidate.alphas, \(al) sum(al / (seq_len(p) - 1 + al)))
      alpha <- candidate.alphas[which.min(abs(prior.K - 3L))]
      rm(candidate.alphas, prior.K); gc()
    } else { alpha <- 1}
  }
  variance.prior.type <- if (is.null(variance.prior.type)) {
    if (prior == "bnp.lasso") "independent" else "conjugate"
  } else { variance.prior.type }
  
  # Run Gibbs sampler ----------------------------------------------------------
  bnplasso_out <- mcmc_sampler(
    X = NULL, y = y, a = a, b = b, alpha = alpha, intercept = FALSE, 
    penalty.type = prior, variance.prior.type = variance.prior.type, 
    max.iters = max.iters, burn.in = burn.in, thin = thin, 
    polya = polya, float = float, sp.means = TRUE 
  )
  
  # Prepare the returns --------------------------------------------------------
  n.obs <- n.preds <- length(y)
  n.draws <- nrow(bnplasso_out$post.beta)
  bnplasso_out[["a"]] <- a
  bnplasso_out[["b"]] <- b
  bnplasso_out[["shrinakge.prior"]] <- prior
  bnplasso_out[["variance.prior.type"]] <- variance.prior.type
  bnplasso_out[["max.iters"]] <- max.iters
  bnplasso_out[["burn.in"]] <- burn.in
  bnplasso_out[["thin"]] <- thin
  bnplasso_out[["n.obs"]] <- n.obs
  bnplasso_out[["n.preds"]] <- n.preds
  bnplasso_out[["n.draws"]] <- n.draws # After burn-in and thinning
  if (prior == "bnp.lasso") bnplasso_out[["alpha"]] <- alpha
  
  # Posterior predictive draws -------------------------------------------------
  post_pred <- loglik <- matrix(data = NA, nrow = n.draws, ncol = n.obs)
  post_sigma <- sqrt(pmax(abs(bnplasso_out$post.sigma2), 1e-16))
  for (s in seq_len(n.draws)) {
    post_beta <- bnplasso_out$post.beta[s,]
    post_pred[s,] <- post_beta + rnorm(n.obs, 0, post_sigma[s])
    loglik[s,] <- dnorm(y, post_beta, post_sigma[s], log = TRUE)
  }
  rm(post_sigma, post_beta); gc()
  bnplasso_out[["y"]] <- y
  bnplasso_out[["loglik"]] <- loglik
  bnplasso_out[["post.pred"]] <- post_pred
  class(bnplasso_out) <- "spmBayes"
  bnplasso_out
}
