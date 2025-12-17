# ------------------------------------------------------------------------------
# S3 methods for class "lmBayes"
# ------------------------------------------------------------------------------

# Class specific helpers -------------------------------------------------------

# Checks if an object inherits from "lmBayes".
is.lmBayes <- \(object) inherits(object, "lmBayes")

# Checks if an object inherits from "spmBayes".
is.spmBayes <- \(object) inherits(object, "spmBayes")

# Summary table for regression coefficients
coef_summary <- function(object, ...) {
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  summary.values <- c("mean", "sd", paste0(probs * 100, "%"), "n_eff")
  beta.names <- paste0("beta_", seq_len(object$n.preds))
  # Summary values for the intercept
  if (!is.spmBayes(object)) {
    if (object$intercept) {
      m <- object$post.mu
      mu.summary <- rep(NA, 8)
      mu.summary[1:7] <- round(c(mean(m), sd(m), quantile(m, probs)), 3)
      mu.summary[8] <- round(coda::effectiveSize(m), 4)
    }
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
  if (!is.spmBayes(object)) {
    if (object$intercept) {
      summary_out <- rbind(mu.summary, summary_out)
      row.names(summary_out) <- c("Intercept", beta.names)
    } else {
      row.names(summary_out) <- beta.names
    }
  } else { row.names(summary_out) <- beta.names }
  colnames(summary_out) <- summary.values
  summary_out
}


# Posterior predictive checks - density plot
postpred.chck.plt <- function(object, ...) {
  n.draws <- object$n.draws
  fitted.vals <- if (!is.spmBayes(object)) {
    object$post.pred.fitted.values
  } else { object$post.pred }
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
  fitted.vals <- if (!is.spmBayes(object)) {
    object$post.pred.fitted.values
  } else { object$post.pred }
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

# Common internal routines for "lmBayes" and "spmBayes" ------------------------

.print_common <- function(x) coef_summary(x)

.coef_common <- function(x) as.data.frame(coef_summary(x)[,-8])

.summary_common <- function(x) {
  if (!is.spmBayes(x)) {
    cat("\n"); cat("\n")
    cat("LINEAR REGRESSION MODEL \n")
  } else {
    cat("\n"); cat("\n")
    cat("SPARSE MEANS PROBLEM \n")
  }
  cat("\n"); cat("\n")
  if (x$shrinakge.prior == "bnp.lasso") {
    cat("NONPARAMETRIC BAYESIAN LASSO \n")
  } else if (x$shrinakge.prior == "b.lasso") {
    cat("BAYESIAN LASSO \n")
  } else if (x$shrinakge.prior == "b.adapt.lasso") {
    cat("BAYESIAN ADAPTIVE LASSO \n")
  }
  cat("\n"); cat("\n")
  cat("Call details: \n")
  cat("\n")
  cat("a =", x$a, "\n")
  cat("b =", x$b, "\n")
  if (x$shrinakge.prior == "bnp.lasso") {
    cat("alpha =", x$alpha, "\n")
  }
  cat("n.obs =", x$n.obs, "\n")
  cat("n.preds =", x$n.preds, "\n")
  cat("n.draws =", x$n.draws, "(after burn-in and thinning) \n")
  cat("elapsed =", format(x$elapsed), "\n")
  cat("\n")
  cat("Coefficients: \n")
  coef_summary(x) 
}

.plot_common <- function(x) {
  # Posterior predictive checks - density plot
  postpred.chck.plt(x)
  readline(prompt = "Press [Enter] to continue...")
  cat()
  # Posterior predictive checks - sample statistics plot
  postpred.samplestats.plt(x)
  # Residual checks plot
  if (!is.spmBayes(x)) {
    readline(prompt = "Press [Enter] to continue...")
    cat()
    res.chck.plt(x)
  }
}

.fitted_common <- function(x) {
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  fv <- if (!is.spmBayes(x)) x$post.pred.fitted.values else x$post.pred
  fv_summary <- apply(fv, 2, \(m) c(mean(m), sd(m), quantile(m, probs))) |> t()
  colnames(fv_summary) <- c("mean", "sd", paste0(probs * 100, "%"))
  as.data.frame(fv_summary)
}

# Class specific main routines -------------------------------------------------


#' Print the results for an object of class \code{lmBayes}
#'
#' @param x An object of class \code{lmBayes}.
#' @param ... Further arguments passed to.
#' 
#' @author Santiago Marin
#'
print.lmBayes <- function(x, ...) {
  if (!is.lmBayes(x)) stop("x should be of class 'lmBayes'")
  .print_common(x)
}

#' Print the results for an object of class \code{spmBayes}
#'
#' @param x An object of class \code{spmBayes}.
#' @param ... Further arguments passed to.
#' 
#' @author Santiago Marin
#'
print.spmBayes <- function(x, ...) {
  if (!is.spmBayes(x)) stop("x should be of class 'spmBayes'")
  .print_common(x)
}


#' Summary table of the results for an object of class \code{lmBayes}
#'
#' @param object An object of class \code{lmBayes}.
#' @param ... Further arguments passed to.
#' 
#' @author Santiago Marin
#'
summary.lmBayes <- function(object, ...) {
  if (!is.lmBayes(object)) stop("object should be of class 'lmBayes'")
  .summary_common(object)
}

#' Summary table of the results for an object of class \code{spmBayes}
#'
#' @param object An object of class \code{spmBayes}.
#' @param ... Further arguments passed to.
#' 
#' @author Santiago Marin
#'
summary.spmBayes <- function(object, ...) {
  if (!is.spmBayes(object)) stop("object should be of class 'spmBayes'")
  .summary_common(object)
}


#' Regression coefficients
#'
#' Extracts the posterior distribution of the regression coefficients from an
#' object of class \code{lmBayes}.
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
  if (!is.lmBayes(object)) stop("object should be of class 'lmBayes'")
  .coef_common(object)
}
coefficients.lmBayes <- coef.lmBayes

#' Mean parameters
#'
#' Extracts the posterior distribution of the mean parameters from an
#' object of class \code{spmBayes}.
#'
#' @aliases coefficients.spmBayes
#'
#' @param object An object of class \code{'spmBayes'}.
#' @param ... Further arguments passed to.
#'
#' @return A data.frame with summary statistics of the posterior distribution of
#' the mean parameters.
#' 
#' @author Santiago Marin
#'
coef.spmBayes <- function(object, ...) {
  if (!is.spmBayes(object)) stop("object should be of class 'spmBayes'")
  .coef_common(object)
}
coefficients.spmBayes <- coef.spmBayes


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
  if (!is.lmBayes(x)) stop("x should be of class 'lmBayes'")
  .plot_common(x)
}

#' Plot the results and diagnostics for an object of class \code{spmBayes}
#'
#' Produces posterior predictive and residual diagnostic plots for an object of 
#' class \code{'spmBayes'}.
#'
#' @param x An object of class \code{'spmBayes'}.
#' @param ... Further arguments passed to.
#' 
#' @author Santiago Marin
#'
plot.spmBayes <- function(x, ...) {
  if (!is.spmBayes(x)) stop("x should be of class 'spmBayes'")
  .plot_common(x)
}


#' Fitted values
#'
#' Extracts the posterior predictive fitted values from an object of class
#' \code{lmBayes}.
#'
#' @aliases fitted.values.lmBayes
#'
#' @param object An object of class \code{'lmBayes'}.
#' @param ... Further arguments passed to.
#'
#' @return A data.frame with summary statistics of the posterior predictive 
#' fitted values.
#' 
#' @author Santiago Marin
#'
fitted.lmBayes <- function(object, ...) {
  if (!is.lmBayes(object)) stop("object should be of class 'lmBayes'")
  .fitted_common(object)
}
fitted.values.lmBayes <- fitted.lmBayes

#' Fitted values
#'
#' Extracts the posterior predictive fitted values from an object of class
#' \code{spmBayes}.
#'
#' @aliases fitted.values.spmBayes
#'
#' @param object An object of class \code{'spmBayes'}.
#' @param ... Further arguments passed to.
#'
#' @return A data.frame with summary statistics of the posterior predictive 
#' fitted values.
#' 
#' @author Santiago Marin
#'
fitted.spmBayes <- function(object, ...) {
  if (!is.spmBayes(object)) stop("object should be of class 'lmBayes'")
  .fitted_common(object)
}
fitted.values.spmBayes <- fitted.spmBayes


#' Residuals from an object of class \code{lmBayes}
#'
#' Extracts the posterior predictive residuals from an object of class
#' \code{lmBayes}.
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
