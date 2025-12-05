# Make the matrix m positive definite ------------------------------------------
make_posdef <- function(m) {
  d <- dim(m)[1]
  eigen_m <- eigen(m, symmetric = TRUE)
  evals_m <- eigen_m$values
  tol <- 2 * (d * max(abs(evals_m)) * .Machine$double.eps)
  delta = pmax(0, tol - evals_m)
  m + eigen_m$vectors %*% diag(delta, d) %*% t(eigen_m$vectors)
}

# Helpers for input validation -------------------------------------------------

# Check if the input y is not a valid response vector
not.y <- \(y) !is.numeric(y) || !is.vector(y) || any(is.na(y))

# Check if the input is not a positive scalar
not.ps <- \(x) (!is.numeric(x)) || (length(x) != 1) || (x <= 0)

# Check if the input is not (1D) logical
not.logic <- \(x) (!is.logical(x)) || (length(x) != 1)

# Check if the input is not a valid variance prior type
not.var <- \(x) !(x %in% c("independent", "conjugate"))

# Check if the input is not a valid shrinkage prior 
not.pr <- \(x) !(x %in% c("bnp.lasso", "b.lasso", "b.adapt.lasso"))

# Check if the input is not (1D) integer
not.int <- \(x) (x %% 1) != 0 || (x <= 0) || (length(x) != 1)

# Check if the a numeric matrix
not.mat <- \(x) !is.numeric(x) || !is.matrix(x) || any(is.na(x))

# Check if the input X is not a valid matrix of covariates
not.x <- \(x, y) {
  !is.numeric(x) || !is.matrix(x) || !(length(y) == nrow(x)) || any(is.na(x))
}
