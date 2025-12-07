# ------------------------------------------------------------------------------
# S3 methods for class "lmBayes"
# ------------------------------------------------------------------------------

# Class specific helpers -------------------------------------------------------

# Checks if an object inherits from "lmBayes".
is.lmBayes <- \(object) inherits(object, "lmBayes")

# Summary table for regression coefficients
coef_summary <- function(object, ...) {
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  summary.values <- c("mean", "sd", paste0(probs * 100, "%"), "n_eff")
  beta.names <- paste0("beta_", seq_len(object$n.preds))
  # Summary values for the intercept
  if (object$intercept) {
    m <- object$post.mu
    mu.summary <- rep(NA, 8)
    mu.summary[1:7] <- round(c(mean(m), sd(m), quantile(m, probs)), 3)
    mu.summary[8] <- round(coda::effectiveSize(m), 4)
  }
  # Summary values for the betas
  summary_out <- apply(
    object$post.beta, 2, \(b) {
      b.summary <- rep(NA, 8)
      b.summary[1:7] <- round(c(mean(b), sd(b), quantile(b, probs)), 3)
      b.summary[8] <- round(coda::effectiveSize(b), 4)
      b.summary
    }
  ) |> t()
  # Prepare the returns
  if (object$intercept) {
    summary_out <- rbind(mu.summary, summary_out)
    row.names(summary_out) <- c("Intercept", beta.names)
  } else {
    row.names(summary_out) <- beta.names
  }
  colnames(summary_out) <- summary.values
  summary_out
}


# Posterior predictive checks - density plot
postpred.chck.plt <- function(object, ...) {
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  n.draws <- object$n.draws
  fitted.vals <- object$post.pred.fitted.values
  x_range <- density(object$y)$x
  x_eval <- seq(x_range[1], x_range[length(x_range)], length.out = 1000)
  y_eval <- matrix(NA, nrow = n.draws, ncol = 1000)
  for (s in seq_len(n.draws)) {
    y_eval[s,] <- FKSUM::fk_density(fitted.vals[s,], x_eval = x_eval)$y
  }
  quants <- apply(y_eval, 2, \(m) quantile(m, c(0.025, 0.975)))
  y.plot_lower <- quants[1,]
  y.plot_upper <- quants[2,]
  par(mfrow = c(1, 1))
  light.col <- "#ADD8E6"
  dark.col <- "#0000FF"
  plot(
    x = x_eval, y = y.plot_upper, type = "n", xlab = "y",
    ylab = "Density", main = "Posterior predictive checks"
  )
  polygon(
    x = c(x_eval, rev(x_eval)), y = c(y.plot_lower, rev(y.plot_upper)),
    col = light.col, border = NA
  )
  lines(
    x = x_eval, y = FKSUM::fk_density(object$y, x_eval = x_eval)$y, 
    col = dark.col, lwd = 2, lty = 2
  )
  legend(
    "topleft", legend = c("y", expression(y[rep])), lty = c(2, 1),
    col = c(dark.col, light.col), lwd = c(2, 15), bty = "n"
  )
  rm(x_range, fitted.vals, quants, y_eval)
  rm(y.plot_lower, y.plot_upper, x_eval); gc()
  par(mfrow = c(1, 1))
}


# Posterior predictive checks - sample statistics plot
postpred.samplestats.plt <- function(object, ...) {
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  fitted.vals <- object$post.pred.fitted.values
  postpred.stats <- apply(fitted.vals, 1, \(m) c(mean(m), sd(m)))
  observed.stats <- c(mean(object$y), sd(object$y))
  pp <- "Posterior predictive "
  mains <- paste0(pp, c("mean", "std. dev."))
  xlabs <- c("mean(y)", "sd(y)")
  par(mfrow = c(1, 2))
  light.col <- "#ADD8E6"
  dark.col <- "#0000FF"
  # Posterior predictive mean - histogram
  postpred.mean.range <- range(postpred.stats[1,])
  min.mean <- min(postpred.mean.range, observed.stats[1])
  max.mean <- max(postpred.mean.range, observed.stats[1])
  hist(
    x = postpred.stats[1,], col = light.col, border = 0, main = mains[1],
    xlim = c(min.mean, max.mean), ylab = "", xlab = xlabs[1], freq = FALSE
  )
  abline(v = observed.stats[1], lwd = 2, col = dark.col, lty = 2)
  legend(
    "topleft", legend = c("y", expression(y[rep])), lty = c(2, 1),
    col = c(dark.col, light.col), lwd = c(2, 15), bty = "n"
  )
  box()
  # Posterior predictive standard deviation - histogram
  postpred.sd.range <- range(postpred.stats[2,])
  min.sd <- min(postpred.sd.range, observed.stats[2])
  max.sd <- max(postpred.sd.range, observed.stats[2])
  hist(
    x = postpred.stats[2,], col = light.col, border = 0, main = mains[2],
    xlim = c(min.sd, max.sd), ylab = "", xlab = xlabs[2], freq = FALSE
  )
  abline(v = observed.stats[2], lwd = 2, col = dark.col, lty = 2)
  legend(
    "topleft", legend = c("y", expression(y[rep])), lty = c(2, 1),
    col = c(dark.col, light.col), lwd = c(2, 15), bty = "n"
  )
  box()
  rm(postpred.stats, mains, fitted.vals, min.mean, min.sd, max.mean, max.sd)
  rm(postpred.mean.range, postpred.sd.range, xlabs, pp); gc()
  par(mfrow = c(1, 1))
}


# Residual checks plot
res.chck.plt <- function(object, ...) {
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  med.res <- apply(object$post.pred.residuals, 2, \(m) median(m, na.rm = T))
  med.fv <- apply(object$post.pred.fitted.values, 2, \(m) median(m, na.rm = T))
  # Internally studentized residuals
  stud.res <- med.res / sd(med.res, na.rm = TRUE)
  par(mfrow = c(1, 1))
  # Residuals plot
  plot(
    x = med.fv, y = med.res, pch = 16, main = "Residuals plot",
    xlab = "Fitted values", ylab = "Residuals"
  )
  lines(lowess(med.fv, med.res), col = 2, lwd = 2, lty = 2)
  readline(prompt = "Press [Enter] to continue...")
  # Studentized residuals plot
  plot(
    x = med.fv, y = stud.res, pch = 16, xlab = "Fitted values", 
    ylab = "Studentized residuals", main = "Studentized residuals plot"
  )
  lines(lowess(med.fv, stud.res), col = 2, lwd = 2, lty = 2)
  rm(med.res, med.fv, stud.res) ; gc()
  par(mfrow = c(1, 1))
}


# Posterior predictive distribution for a new data set 'X.new' of 
# dimension n.new \times p, where each row is a test observation.
#
# Returns a matrix, where each row corresponds to an MCMC draw and 
# each column to a test (held-out) observation.
postpred.newdata <- function(object, X.new) {
  n.test <- nrow(X.new)
  n.draws <- object$n.draws
  # Pre-compute all linear predictors (without the intercept)
  linPreds <- tcrossprod(X.new, object$post.beta)
  # Add the intercept if needed
  if (object$intercept) {
    linPreds <- linPreds + matrix(
      data = object$post.mu, nrow = nrow(linPreds),
      ncol = ncol(linPreds), byrow = TRUE
    )
  }
  # Note: In "post_pred_fits" each row corresponds to an
  # MCMC draw and each column to a test (held-out) observation.
  post_pred_fits <- matrix(data = NA, nrow = n.draws, ncol = n.test)
  post_sigmas <- sqrt(abs(object$post.sigma2) + 1e-16)
  for (s in seq_len(n.draws)) {
    fit_val_s <- linPreds[,s] + rnorm(n.test, 0, post_sigmas[s])
    post_pred_fits[s,] <- fit_val_s
  }
  rm(linPreds, post_sigmas, fit_val_s); gc()
  post_pred_fits
}


# Class specific main routines -------------------------------------------------


#' Print the results for an object of class \code{'lmBayes'}
#'
#' @param x An object of class \code{'lmBayes'}.
#' @param ... Further arguments passed to.
#' 
#' @author Santiago Marin
#'
print.lmBayes <- function(x, ...) {
  chck <- !is.lmBayes(x)
  if (chck) stop("object should be of class 'lmBayes'")
  coef_summary(x)
}


#' Summary table of the results for an object of class \code{'lmBayes'}
#'
#' @param object An object of class \code{'lmBayes'}.
#' @param ... Further arguments passed to.
#' 
#' @author Santiago Marin
#'
summary.lmBayes <- function(object, ...) {
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  cat("\n"); cat("\n")
  if (object$shrinakge.prior == "bnp.lasso") {
    cat("NONPARAMETRIC BAYESIAN LASSO \n")
  } else if (object$shrinakge.prior == "b.lasso") {
    cat("BAYESIAN LASSO \n")
  } else if (object$shrinakge.prior == "b.adapt.lasso") {
    cat("BAYESIAN ADAPTIVE LASSO \n")
  }
  cat("\n"); cat("\n")
  cat("Call details: \n")
  cat("\n")
  cat("a =", object$a, "\n")
  cat("b =", object$b, "\n")
  if (object$shrinakge.prior == "bnp.lasso") {
    cat("alpha =", object$alpha, "\n")
  }
  cat("n.obs =", object$n.obs, "\n")
  cat("n.preds =", object$n.preds, "\n")
  cat("n.draws =", object$n.draws, "(after burn-in and thinning) \n")
  cat("elapsed =", format(object$elapsed), "\n")
  cat("\n")
  cat("Coefficients: \n")
  coef_summary(object) 
}


#' Plot the results and diagnostics for an object of class \code{lmBayes}
#'
#' Produces posterior predictive and residual diagnostic plots for an object of 
#' class \code{'lmBayes'}.
#'
#' @param x An object of class \code{'lmBayes'}.
#' @param ... Further arguments passed to.
#' 
#' @author Santiago Marin
#'
plot.lmBayes <- function(x, ...) {
  chck <- !is.lmBayes(x)
  if (chck) stop("object should be of class 'lmBayes'")
  # Posterior predictive checks - density plot
  postpred.chck.plt(x)
  readline(prompt = "Press [Enter] to continue...")
  cat()
  # Posterior predictive checks - sample statistics plot
  postpred.samplestats.plt(x)
  readline(prompt = "Press [Enter] to continue...")
  cat()
  # Residual checks plot
  res.chck.plt(x)
}


#' Fitted Values from an object of class \code{lmBayes} 
#'
#' Summary statistics of the posterior predictive distribution of the 
#' fitted values from an object of class \code{lmBayes}.
#'
#' @aliases fitted.values.lmBayes
#'
#' @param object An object of class \code{'lmBayes'}.
#' @param ... Further arguments passed to.
#'
#' @return A data.frame with summary statistics of the posterior predictive 
#' distribution of the fitted values.
#' 
#' @author Santiago Marin
#'
fitted.lmBayes <- function(object, ...) {
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  fitted.summary <- apply(object$post.pred.fitted.values, 2, \(m) {
    c(mean(m), sd(m), quantile(m, probs))
  }) |> t()
  colnames(fitted.summary) <- c("mean", "sd", paste0(probs * 100, "%"))
  as.data.frame(fitted.summary)
}
fitted.values.lmBayes <- fitted.lmBayes


#' Residuals from an object of class \code{lmBayes}
#'
#' Summary statistics of the posterior predictive distribution of the 
#' residuals from an object of class \code{lmBayes}.
#'
#' @aliases resid.lmBayes
#'
#' @param object An object of class \code{'lmBayes'}.
#' @param ... Further arguments passed to.
#'
#' @return A data.frame with summary statistics of the posterior predictive 
#' distribution of the model residuals.
#' 
#' @author Santiago Marin
#'
residuals.lmBayes <- function(object, ...) {
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  res.summary <- apply(object$post.pred.residuals, 2, \(m) {
    c(mean(m), sd(m), quantile(m, probs))
  }) |> t()
  colnames(res.summary) <- c("mean", "sd", paste0(probs * 100, "%"))
  as.data.frame(res.summary)
}
resid.lmBayes <- residuals.lmBayes


#' Regression coefficients from an object of class \code{lmBayes}
#'
#' Summary statistics of the posterior distribution of the regression 
#' coefficients from an object of class \code{lmBayes}.
#'
#' @aliases coefficients.lmBayes
#'
#' @param object An object of class \code{'lmBayes'}.
#' @param ... Further arguments passed to.
#'
#' @return A data.frame with summary statistics of the posterior distribution of
#' the regression coefficients.
#' 
#' @author Santiago Marin
#'
coef.lmBayes <- function(object, ...) {
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  out <- coef_summary(object)
  out <- out[,-8]
  as.data.frame(out)
}
coefficients.lmBayes <- coef.lmBayes


#' Posterior predictive distribution for new data
#'
#' Compute the posterior predictive distribution for an object of class 
#' \code{lmBayes}.
#'
#' @param object An object of class \code{'lmBayes'}.
#' @param X.new A new matrix of predictors, where each row is a new observation.
#' @param ... Further arguments passed to. 
#'
#' @return A matrix where each row corresponds to an MCMC draw and each column 
#' to an observation in the new data.
#'
#' @author Santiago Marin
#'
predict.lmBayes <- function(object, X.new, ...) {
  chck <- !is.lmBayes(object)
  if (chck) stop("object should be of class 'lmBayes'")
  n.draws <- object$n.draws
  n.preds <- object$n.preds
  if(not.mat(X.new)) stop("X.new must be a valid numeric matrix")
  if(ncol(X.new) != n.preds) stop("X.new has an incorrect number of columns")
  postpred.newdata(object, X.new)
}
