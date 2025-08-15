# Make the matrix m positive definite ------------------------------------------
make_posdef <- function(m) {
  d <- dim(m)[1]
  eigen_m <- eigen(m, symmetric = TRUE)
  evals_m <- eigen_m$values
  tol <- 2 * (d * max(abs(evals_m)) * .Machine$double.eps)
  delta = pmax(0, tol - evals_m)
  m + eigen_m$vectors %*% diag(delta, d) %*% t(eigen_m$vectors)
}
