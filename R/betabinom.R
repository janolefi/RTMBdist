#' Beta-binomial distribution
#'
#' Density and random generation for the beta-binomial distribution.
#'
#' @details
#' This implementation of \code{dbetabinom} allows for automatic differentiation with \code{RTMB}.
#'
#' \deqn{P(X = k;\, n, a, b) = \binom{n}{k} \frac{B(k+a,\, n-k+b)}{B(a,\, b)}, \quad k = 0, 1, \ldots, n.}
#'
#' @param x vector of non-negative counts.
#' @param size vector of total counts (number of trials). Needs to be >= \code{x}.
#' @param n number of random values to return (for \code{rbetabinom}).
#' @param shape1 positive shape parameter 1 of the Beta prior.
#' @param shape2 positive shape parameter 2 of the Beta prior.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#'
#' @return
#' \code{dbetabinom} gives the density and \code{rbetabinom} generates random samples.
#'
#' @examples
#' set.seed(123)
#' x <- rbetabinom(1, 10, 2, 5)
#' d <- dbetabinom(x, 10, 2, 5)
#' @name betabinom
NULL

#' @rdname betabinom
#' @export
#' @import RTMB
dbetabinom <- function(x, size, shape1, shape2, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # if (any(x < 0) || any(x != floor(x)))
    #   stop("x must be non-negative integers.")
    if (any(size < 0) || any(size != floor(size)))
      stop("size must be non-negative integers.")
    # if (any(x > size))
    #   stop("x cannot be greater than size.")
    if (any(shape1 <= 0) || any(shape2 <= 0))
      stop("shape1 and shape2 must be positive.")
  }

  if (inherits(x, "simref")) {
    return(dGenericSim("dbetabinom", x = x, size = size, shape1 = shape1, shape2 = shape2, log = log))
  }
  if (inherits(x, "osa")) {
    stop("Beta-binomial does not support OSA residuals.")
  }

  # recycle scalars
  nx <- length(x)
  if (length(size) == 1) size <- rep(size, nx)
  if (length(shape1) == 1) shape1 <- rep(shape1, nx)
  if (length(shape2) == 1) shape2 <- rep(shape2, nx)

  logdens <- lgamma(size + 1) - lgamma(x + 1) - lgamma(size - x + 1) +
    lgamma(x + shape1) + lgamma(size - x + shape2) -
    lgamma(size + shape1 + shape2) +
    lgamma(shape1 + shape2) - lgamma(shape1) - lgamma(shape2)

  if (log) return(logdens)
  exp(logdens)
}
#' @rdname betabinom
#' @export
rbetabinom <- function(n, size, shape1, shape2) {

  if (length(size) == 1) size <- rep(size, n)
  if (length(shape1) == 1) shape1 <- rep(shape1, n)
  if (length(shape2) == 1) shape2 <- rep(shape2, n)

  if (any(shape1 <= 0) || any(shape2 <= 0))
    stop("shape1 and shape2 must be positive.")

  # sample success probabilities from Beta(shape1, shape2)
  p <- rbeta(n, shape1, shape2)

  # sample from Binomial(size, p)
  stats::rbinom(n, size, p)
}
