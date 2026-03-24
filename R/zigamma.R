#' Zero-inflated gamma distribution
#'
#' Density, distribution function, and random generation for
#' the zero-inflated gamma distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return
#' @param shape positive shape parameter
#' @param scale positive scale parameter
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzigamma} gives the density, \code{pzigamma} gives the distribution function, and \code{rzigamma} generates random deviates.
#'
#' @examples
#' x <- rzigamma(1, 1, 1, 0.5)
#' d <- dzigamma(x, 1, 1, 0.5)
#' p <- pzigamma(x, 1, 1, 0.5)
#' @name zigamma
NULL

#' @rdname zigamma
#' @export
#' @importFrom RTMB dgamma logspace_add
dzigamma = function(x, shape, scale, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure shape >= 0, scale > 0, zeroprob in [0,1]
    if (any(shape < 0)) stop("shape must be >= 0")
    if (any(scale <= 0)) stop("scale must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dzigamma", x=x, shape = shape, scale = scale, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzigamma", x=x, shape = shape, scale = scale, zeroprob = zeroprob, log=log))
  }

  eps <- .Machine$double.xmin # so that gradient is not NaN bc -Inf * 0

  logdens <- RTMB::dgamma(x + eps, shape = shape, scale = scale, log = TRUE)
  logdens <- log_zi(x, logdens, zeroprob)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname zigamma
#' @importFrom RTMB pgamma
#' @export
pzigamma <- function(q, shape, scale, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure shape >= 0, scale > 0, zeroprob in [0,1]
    if (any(shape < 0)) stop("shape must be >= 0")
    if (any(scale <= 0)) stop("scale must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # s1 <- 2 * sign(q) - 1 # gives -3 for q < 0, -1 for q == 0, and 1 for q > 0
  # s2 <- sign(3 + s1) # only zero or 1

  # cdf <- 0.5 * (1 - s1) * zeroprob +
  #   0.5 * (1 + s1) * (zeroprob + (1 - zeroprob) * pgamma(q, shape, scale))
  # cdf <- cdf * s2 # set negative values to 0

  p <- iszero(q) * zeroprob +
    ispos_strict(q) * (zeroprob + (1 - zeroprob) * pgamma(q, shape, scale))

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname zigamma
#' @importFrom stats runif rgamma
#' @export
rzigamma <- function(n, shape, scale, zeroprob = 0) {
  # ensure shape >= 0, scale > 0, zeroprob in [0,1]
  if (any(shape < 0)) stop("shape must be >= 0")
  if (any(scale <= 0)) stop("scale must be > 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  u <- runif(n)
  res <- rep(0, n)
  is_zero <- u < zeroprob
  res[!is_zero] <- rgamma(sum(!is_zero), shape = shape, scale = scale)

  return(res)
}
