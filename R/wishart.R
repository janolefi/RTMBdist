#' Wishart distribution
#'
#' Density and random generation for the wishart distribution
#'
#' @param x positive definite \eqn{p \times p} matrix or array of such matrices of dimension \eqn{p \times p \times n} (for \eqn{n} density evaluations)
#' @param nu degrees of freedom, needs to be greater than \code{p - 1}
#' @param Sigma scale matrix, needs to be positive definite and match the dimension of \code{x}.
#' @param log logical; if \code{TRUE}, densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param n number of random deviates to return
#'
#' @returns
#' \code{dwishart} gives the density,
#' \code{rwishart} generates random deviates (matrix for \code{n = 1}, array with \code{n} slices for \code{n > 1})
#'
#' @examples
#' # single input: matrix-valued
#' x <- rwishart(1, nu = 5, Sigma = diag(3))
#' d <- dwishart(x, nu = 5, Sigma = diag(3))
#'
#' # multiple inputs: array of matrices
#' x <- rwishart(4, nu = 5, Sigma = diag(3))
#' d <- dwishart(x, nu = 5, Sigma = diag(3))
#'
#' # multiple inputs for x, nu and Sigma
#' nu <- c(7,5,8,9)
#' Sigma <- array(dim = c(3,3,4))
#' for(i in 1:4) Sigma[,,i] <- (i + 3) * diag(3)
#' x <- rwishart(4, nu, Sigma)
#' d <- dwishart(x, nu, Sigma)
#' @name wishart
NULL

#' @rdname wishart
#' @export
#' @import RTMB
dwishart <- function(x, nu, Sigma, log = FALSE) {

  if(inherits(x, "simref")) {
    if(is.matrix(x)) {
      n <- 1 # if x is matrix, only one single matrix-valued sample
    } else if(is.array(x)) {
      n <- dim(x)[3] # if x is array, samples are in third dimension
    }
    x[] <- rwishart(n, nu=nu, Sigma=Sigma)
    return(0)
  }
  if(inherits(x, "osa")) {
    stop("Wishart does not support OSA residuals.")
  }

  # Handle array input ------------------------------------------------------

  if(length(dim(x)) == 3) {

    p <- dim(x)[1] # matrix dimension
    n <- dim(x)[3] # number of evaluation points

    if(dim(x)[2] != p) { # square check
      stop("x must contain square matrices")
    }

    # broadcast nu
    if(length(nu) == 1) nu <- rep(nu, n)
    if(length(nu) != n) stop("nu must have length 1 or n")

    # broadcast Sigma
    if(length(dim(Sigma)) == 2) {
      Sigma <- AD(array(Sigma, dim = c(p, p, n))) # AD wrap
    }
    if(!all(dim(Sigma) == c(p,p,n))) {
      stop("Sigma must be 'p x p' or 'p x p x n'")
    }

    res <- sapply(seq_len(n),
                  function(i) dwishart(x[,,i], nu[i], Sigma[,,i], log = TRUE))

    if(log) return(res)
    return(exp(res))
  }


  # Single matrix case ------------------------------------------------------

  p <- nrow(x) # matrix dimension

  if(ncol(x) != p) { # square check
    stop("x must contain square matrices")
  }

  # Sigma has correct dim?
  if (!all(dim(x) == dim(Sigma))) {
    stop("x and Sigma must have the same dimensions")
  }

  # parameter-dependent checks only in non-AD mode
  if(!ad_context()) {

    args <- as.list(environment())
    simulation_check(args)

    if (!isSymmetric(x)) stop("x must be symmetric")
    if (!isSymmetric(Sigma)) stop("Sigma must be symmetric")
    if (nu <= p - 1) stop("nu must be greater than p - 1")
  }

  logdet_Sigma <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  logdet_x <- as.numeric(determinant(x, logarithm = TRUE)$modulus)

  logC <- - (nu * p / 2) * log(2) - (nu / 2) * logdet_Sigma - lmultigamma(nu/2, p)

  logdens <- logC +
    ((nu - p - 1)/2) * logdet_x -
    0.5 * sum(solve(Sigma) * x)

  if(log) return(logdens)
  exp(logdens)
}


#' @rdname wishart
#' @export
#' @importFrom stats rchisq rnorm
rwishart <- function(n, nu, Sigma) {

  # determine dimension
  if(length(dim(Sigma)) == 2) {
    p <- nrow(Sigma)
    if(ncol(Sigma) != p) stop("Sigma must be square")
  } else if(length(dim(Sigma)) == 3) {
    p <- dim(Sigma)[1]
    if(dim(Sigma)[2] != p) stop("Sigma must contain square matrices")
  } else {
    stop("Sigma must be 'p x p' or 'p x p x n'")
  }

  # broadcast nu
  if(length(nu) == 1) nu <- rep(nu, n)
  if(length(nu) != n) stop("nu must have length 1 or n")

  # broadcast Sigma
  if(length(dim(Sigma)) == 2) {
    Sigma <- array(Sigma, dim = c(p, p, n))
  }

  if(!all(dim(Sigma) == c(p,p,n))) {
    stop("Sigma must be 'p x p' or 'p x p x n'")
  }

  # generate samples
  res <- lapply(seq_len(n), function(i) rwishart1(nu[i], Sigma[,,i]))

  if (n == 1) return(res[[1]])
  simplify2array(res)
}

# internal for n = 1 - not exported
rwishart1 <- function(nu, Sigma) {

  p <- nrow(Sigma)

  # checks
  if(ncol(Sigma) != p) {
    stop("Sigma must be square")
  }
  if (!isSymmetric(Sigma)) {
    stop("Sigma must be symmetric")
  }
  if (nu <= p - 1) {
    stop("nu must be greater than p - 1")
  }
  if (any(eigen(Sigma, only.values = TRUE)$values <= 0)) {
    stop("Sigma must be positive definite")
  }

  L <- chol(Sigma) # upper-triangular Cholesky
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
  (res + t(res)) / 2
}
