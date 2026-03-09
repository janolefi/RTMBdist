#' Wishart distribution
#'
#' Density and random generation for the wishart distribution
#'
#' @param x positive definite p x p matrix of evaluation points
#' @param nu degrees of freedom, needs to be greater than \code{p - 1}
#' @param Sigma scale matrix, needs to be positive definite and match the dimension of \code{x}.
#' @param log logical; if \code{TRUE}, densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param n number of random deviates to return
#'
#' @returns
#' \code{dwishart} gives the density,
#' \code{rwishart} generates random deviates (matrix for \code{n = 1}, array for \code{n > 1})
#'
#' @examples
#' x <- rwishart(1, nu = 5, Sigma = diag(3))
#' d <- dwishart(x, nu = 5, Sigma = diag(3))
#' @name wishart
NULL

#' @rdname wishart
#' @export
#' @import RTMB
dwishart <- function(x, nu, Sigma, log = FALSE) {

  # potentially escape to RNG or produce error for CDF
    if(inherits(x, "simref")) {
      # n <- if (is.matrix(x)) nrow(x) else 1
      # x[] <- rwishart(n, nu=nu, Sigma=Sigma)
      # return(0)
      stop("Currently, Wishart does not support automatic simulation.")
    }
    if(inherits(x, "osa")) {
      stop("Wishart does not support OSA residuals.")
    }

  p <- nrow(x)

  # Input dimension check
  if (!all(dim(x) == dim(Sigma))) stop("x and Sigma must have the same dimensions")

  # Input checks when not in AD context
  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order

    # Check that x is symmetric
    if (!isSymmetric(x)) stop("x must be a symmetric matrix")
    # Check that Sigma is symmetric
    if (!isSymmetric(Sigma)) stop("Sigma must be a symmetric matrix")
    # Check that nu is valid
    if (nu <= p - 1) stop("nu must be greater than p - 1")
    # Check that Sigma is positive definite
    # if (any(eigen(Sigma, only.values = TRUE)$values <= 0)) stop("Sigma must be positive definite")
    # Check that x is positive definite
    # if (any(eigen(x, only.values = TRUE)$values <= 0)) stop("x must be positive definite")
  }

  # Cholesky decomposition for log-determinant
  logdet_Sigma <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  logdet_x <- as.numeric(determinant(x, logarithm = TRUE)$modulus)

  # Normalizing constant
  logC <- - (nu * p / 2) * log(2) - (nu / 2) * logdet_Sigma - lmultigamma(nu/2, p)

  # Log-likelihood
  logdens <- logC + ((nu - p - 1)/2) * logdet_x - 0.5 * sum(solve(Sigma) * x)

  if(log) return(logdens)
  return(exp(logdens))
}
# dwishart <- function(x, nu, Sigma, log = FALSE) {
#   # potentially escape to RNG or produce error for CDF
#   if(inherits(x, "simref")) {
#     return(dGenericSim("dwishart", x=x, nu=nu, Sigma=Sigma, log=log))
#   }
#   if(inherits(x, "osa")) {
#     stop("Wishart does not support OSA residuals.")
#   }
#
#   if(is.matrix(x)) {
#     # single evaluation matrix
#     if(length(nu) > 1) stop("For a single matrix x, nu can only have length 1.")
#     if(length(dim(Sigma)) != 2) stop("For a single matrix x, Sigma must be a matrix.")
#     return(dwishart1(x, nu, Sigma, log = log))
#   } else if(is.array(x)) {
#     # evaluation on array
#     n <- dim(x)[3]
#     if(length(nu) == 1) nu <- rep(nu, n)
#     if(length(nu) != n) stop("nu must be of length 1 or match the number of samples in x.")
#     if(is.matrix(Sigma)) Sigma <- array(Sigma, dim = c(dim(Sigma), n))
#     if(is.array(Sigma)) {
#       if(dim(Sigma)[3] != n) stop("If Sigma is an array, its number of slices must be 1 or match x.")
#     }
#     res <- sapply(seq_len(n),
#                   function(i) dwishart1(x[,,i], nu[i], Sigma[,,i], log = log))
#   }
# }

#' @rdname wishart
#' @export
#' @importFrom stats rchisq rnorm
rwishart <- function(n, nu, Sigma) {
  p <- nrow(Sigma)
  if (!all(dim(Sigma) == c(p, p))) stop("Sigma must be square")

  # Input checks
  if (!isSymmetric(Sigma)) stop("Sigma must be symmetric")
  if (nu <= p - 1) stop("nu must be greater than p - 1")
  if (any(eigen(Sigma, only.values = TRUE)$values <= 0)) stop("Sigma must be positive definite")

  res <- lapply(seq_len(n), function(i) rwishart1(nu, Sigma))

  if (n == 1) return(res[[1]])
  simplify2array(res)
}

# internal for n = 1 - not exported
rwishart1 <- function(nu, Sigma) {
  p <- nrow(Sigma)

  L <- chol(Sigma)  # upper-triangular Cholesky
  A <- matrix(0, p, p)

  for (i in seq_len(p)) {
    # Diagonal: sqrt of chi-squared
    A[i, i] <- sqrt(stats::rchisq(1, df = nu - i + 1))
    if (i < p) {
      # Lower-triangular off-diagonal: standard normal
      A[(i + 1):p, i] <- stats::rnorm(p - i)
    }
  }

  # Bartlett decomposition: X = L %*% A %*% t(A) %*% t(L)
  res <- L %*% A %*% t(A) %*% t(L)

  # force symmetric
  res <- (res + t(res)) / 2

  return(res)
}


