#' Dirichlet-multinomial distribution
#'
#' Density and and random generation for the Dirichlet-multinomial distribution.
#'
#' @details
#' This implementation of \code{ddirmult} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x vector or matrix of non-negative counts, where rows are observations and columns are categories.
#' @param size vector of total counts for each observation. Needs to match the row sums of \code{x}.
#' @param n number of random values to return.
#' @param alpha vector or matrix of positive shape parameters
#' @param log logical; if \code{TRUE}, densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{ddirmult} gives the density and \code{rdirmult} generates random samples.
#'
#' @examples
#' # single alpha
#' alpha <- c(1,2,3)
#' size <- 10
#' x <- rdirmult(1, size, alpha)
#' d <- ddirmult(x, size, alpha)
#' # vectorised over alpha and size
#' alpha <- rbind(alpha, 2*alpha)
#' size <- c(size, 3*size)
#' x <- rdirmult(2, size, alpha)
#  d <- ddirmult(x, size, alpha)
#' @name dirmult
NULL

#' @rdname dirmult
#' @export
#' @import RTMB
ddirmult <- function(x, size, alpha, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
  }

  # potentially escape to RNG or produce error for CDF
  if(inherits(x, "simref")) {
    n <- if (is.matrix(x)) nrow(x) else 1
    x[] <- rdirmult(n, size=size, alpha=alpha)
    return(0)
  }
  if(inherits(x, "osa")) {
    stop("Dirichlet-multinomial does not support OSA residuals.")
  }

  # Check if x and alpha are vectors by checking if they have dimensions
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  if (is.null(dim(alpha))) alpha <- matrix(alpha, nrow = 1)

  nx <- nrow(x) # number of evaluation points

  # Ensure x and alpha have the same number of columns
  if (ncol(x) != ncol(alpha)) {
    stop("x and alpha must have the same number of columns (categories).")
  }
  # If x is a matrix, check that alpha either has the same number of rows or one
  if (nx > 1) {
    if (nrow(alpha) == 1) {
      alpha <- alpha[rep(1, nx), , drop = FALSE]  # explicit recycle
    } else {
      if (nrow(alpha) != nx) {
        stop("If x is a matrix, alpha must have either one row or the same number of rows as x.")
      }
    }
  }
  # Check if n is a single value or a vector of length
  if (length(size) == 1) {
    size <- rep(size, nx)
  } else if (length(size) != nx) {
    stop("size must be a single value or a vector of length equal to the number of rows in x.")
  }

  alpha0 <- rowSums(alpha)

  ## parameter-dependent checks if not in adcontext
  if(!ad_context()) {
    # Check that all elements of x are non-negative integers
    if (any(x < 0) || any(x != floor(x))) {
      stop("All elements of x must be non-negative integers.")
    }
    # Check that counts sum to size
    if (any(rowSums(x) != size)) {
      stop("Each row of x must sum to the corresponding size.")
    }
    # Check that all elements of alpha are positive
    if (any(alpha <= 0)) {
      stop("All elements of alpha must be positive.")
    }
  }

  logconst <- lgamma(alpha0) + lgamma(size + 1) - lgamma(size + alpha0)

  logdens <- rowSums(lgamma(x + alpha) - lgamma(alpha) - lgamma(x + 1)) +
    logconst

  logdens <- as.numeric(logdens)

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname dirmult
#' @export
#' @importFrom stats rmultinom
rdirmult <- function(n, size, alpha) {

  # Check if alpha is vector by checking if it has dimensions
  if (is.vector(alpha)) alpha <- matrix(alpha, nrow = 1)

  if (n > 1) {
    if (nrow(alpha) == 1) {
      alpha <- alpha[rep(1, n), , drop = FALSE]  # explicit recycle
    } else {
      if (nrow(alpha) != n) {
        stop("alpha must have either 1 row or n rows")
      }
    }
  }
  # Check if n is a single value or a vector of length
  if (length(size) == 1) {
    size <- rep(size, n)
  } else if (length(size) != n) {
    stop("size must be a single value or a vector of length equal to n.")
  }

  # check that all elements of alpha are positive
  if (any(alpha <= 0)) {
    stop("All elements of alpha must be positive.")
  }

  # Draw n sets of theta from Dirichlet(alpha)
  theta <- rdirichlet(n, alpha)

  # Draw multinomial counts
  counts <- t(vapply(seq_len(n), function(i)
    as.vector(stats::rmultinom(1, size[i], theta[i, ])), numeric(ncol(alpha))))

  # colnames(counts) <- paste0("cat", 1:ncol(counts))
  counts
}
