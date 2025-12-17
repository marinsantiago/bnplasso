# ------------------------------------------------------------------------------
# Compute point estimates of the regression coefficients (and the intercept)
# ------------------------------------------------------------------------------

# Point estimate helpers -------------------------------------------------------

# Compute the posterior mean for all the regression coefficients
post.mean <- function(object) colMeans(object$post.beta)

# Compute the posterior mode for a sequence of posterior draws
single.post.mode <- function(x) {
  range_x <- range(x)
  x_eval <- seq(range_x[1], range_x[2], length.out = 2000)
  dens.eval <- FKSUM::fk_density(x, x_eval = x_eval)$y
  x_eval[which.max(dens.eval)]
}

# Compute the posterior mode for all the regression coefficients
post.mode <- function(object) {
  n.preds <- object$n.preds
  sapply(seq_len(n.preds), \(j) single.post.mode(object$post.beta[,j]))
}

# If the 50% credible interval contains zero, remove the predictor. Otherwise,
# retain either the posterior mode or the posterior mean.
credinterval.criterion <- function(object, retain) {
  n.preds <- object$n.preds
  # .5 credible sets
  alpha_half <- 0.5 / 2
  probs <- c(alpha_half, 1 - alpha_half)
  quant <- apply(object$post.beta, 2, \(c) quantile(c, probs = probs))
  lower <- quant[1,]
  upper <- quant[2,]
  sapply(
    seq_len(n.preds), \(j) {
      # Check if the cred. interval contains zero
      out <- if ((lower[j] < 0) && (upper[j] > 0)) {
        0L # Exclude the variable
      } else {
        if (retain == "mean") {
          mean(object$post.beta[,j]) # Retain the posterior mean
        } else if (retain == "mode") {
          single.post.mode(object$post.beta[,j]) # Retain the posterior mode
        }
      }
    }
  )
}

# If the posterior probability of beta_j being in the scaled neighborhood 
# from Li and Lin (2010) is larger than 0.5, remove the predictor. Otherwise, 
# retain either the posterior mode or the posterior mean.
scaled.neighbor.criterion <- function(object, retain) {
  n.preds <- object$n.preds
  n.draws <- object$n.draws
  sapply(
    seq_len(n.preds), \(j) {
      beta_j <- object$post.beta[,j]
      sd.beta_j <- sd(beta_j)
      bounds.beta_j <- c(-sd.beta_j, sd.beta_j)
      interval.chck <- (beta_j > bounds.beta_j[1]) & (beta_j < bounds.beta_j[2])
      prob.beta_j <- sum(interval.chck) / n.draws
      if (prob.beta_j >= 0.5) {
        0L
      } else {
        if (retain == "mean") {
          mean(object$post.beta[,j]) # Retain the posterior mean
        } else if (retain == "mode") {
          single.post.mode(object$post.beta[,j]) # Retain the posterior mode
        }
      }
    }
  )
}


# Point estimates main routine -------------------------------------------------


#' Point estimates of the regression coefficients
#'
#' Computes point estimates of the regression coefficients for an object of 
#' class \code{lmBayes} or \code{'spmBayes'}.
#'
#' @param object An object of class \code{'lmBayes'} or \code{'spmBayes'}.
#' @param type A character string denoting which algorithm should be used to
#'   recover the point estimates. The options are: (i) \code{"sn"} for the
#'   scaled neighborhood criterion from Li and Lin (2010), (ii) 
#'   \code{"cred.int"} which excludes a predictor if the 50\% credible interval 
#'   contains zero, (iii) \code{"post.mode"} for the posterior mode, or (iv)
#'   \code{"post.mean"} for the posterior mean. Default is \code{"sn"}.
#' @param retain A character string denoting which posterior summary should be 
#'   retained if \code{type = "sn"} or \code{type = "cred.int"}. The options 
#'   are: (i) \code{"mode"} to retain the posterior mode or (ii) \code{"mean"} 
#'   to retain the posterior mean.
#'
#' @return A numeric vector containing the point estimates of the regression
#' coefficients.
#' 
#' @references
#'
#' Q. Li, and N. Lin (2010), The Bayesian elastic net. \emph{Bayesian Analysis},
#' 5(1):151-170.
#'
#' @author Santiago Marin
#'
point.estimates <- function(object, type = "sn", retain = "mean") {
  # Input validation -----------------------------------------------------------
  if (!is.lmBayes(object) & !is.spmBayes(object)) {
    stop("object should be of class 'lmBayes' or 'spmBayes'")
  }
  chck <- !(type %in% c("sn", "cred.int", "post.mode", "post.mean"))
  mssg <- "type must be either 'sn', 'cred.int', 'post.mode', or 'post.mean'"
  if (chck) stop(mssg)
  chck <- !(retain %in% c("mode", "mean"))
  if (chck) stop("retain must be either 'mode' or 'mean'")
  
  # Point estimates for beta ---------------------------------------------------
  beta.hat <- if (type == "sn") {
    scaled.neighbor.criterion(object, retain)
  } else if (type == "cred.int") {
    credinterval.criterion(object, retain)
  } else if (type == "post.mode") {
    post.mode(object)
  } else if (type == "post.mean") {
    post.mean(object)
  }
  names(beta.hat) <- paste0("beta_", seq_len(object$n.preds))

  # Point estimates for the intercept ------------------------------------------
  if (!is.spmBayes(object)) {
    if (object$intercept) {
      mu.hat <- if (type == "post.mode") {
        single.post.mode(object$post.mu)
      } else if (type == "post.mean") {
        mean(object$post.mu)
      } else if (retain == "mode") {
        single.post.mode(object$post.mu)
      } else if (retain == "mean") {
        mean(object$post.mu)
      }
      beta.hat <- c(mu.hat, beta.hat)
      names(beta.hat)[1] <- "Intercept"
    }
  }
  beta.hat
}
