#' Zero-inflated Weibull distribution
#'
#' Density, distribution function, and random generation for
#' the zero-inflated Weibull distribution.
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
#' \code{dziweibull} gives the density, \code{pziweibull} gives the distribution function, and \code{rziweibull} generates random deviates.
#'
#' @examples
#' x <- rziweibull(1, 1, 1, 0.5)
#' d <- dziweibull(x, 1, 1, 0.5)
#' p <- pziweibull(x, 1, 1, 0.5)
#' @name ziweibull
NULL

#' @rdname ziweibull
#' @export
#' @importFrom RTMB dweibull logspace_add
dziweibull <- function(x, shape, scale, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure shape > 0, scale > 0, zeroprob in [0,1]
    if (any(shape <= 0)) stop("shape must be > 0")
    if (any(scale <= 0)) stop("scale must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dziweibull", x=x, shape = shape, scale = scale, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dziweibull", x=x, shape = shape, scale = scale, zeroprob = zeroprob, log=log))
  }

  eps <- .Machine$double.xmin # so that gradient is not NaN bc -Inf * 0
  logdens <- RTMB::dweibull(x + eps, shape = shape, scale = scale, log = TRUE)
  logdens <- log_zi(x, logdens, zeroprob)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname ziweibull
#' @importFrom RTMB pweibull
#' @export
pziweibull <- function(q, shape, scale, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure shape > 0, scale > 0, zeroprob in [0,1]
    if (any(shape <= 0)) stop("shape must be > 0")
    if (any(scale <= 0)) stop("scale must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  p <- iszero(q) * zeroprob +
    ispos_strict(q) * (zeroprob + (1 - zeroprob) * pweibull(q, shape, scale))

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname ziweibull
#' @importFrom stats runif rweibull
#' @export
rziweibull <- function(n, shape, scale, zeroprob = 0) {
  # ensure shape > 0, scale > 0, zeroprob in [0,1]
  if (any(shape <= 0)) stop("shape must be > 0")
  if (any(scale <= 0)) stop("scale must be > 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  u <- runif(n)
  res <- rep(0, n)
  is_zero <- u < zeroprob
  res[!is_zero] <- rweibull(sum(!is_zero), shape = shape, scale = scale)

  return(res)
}
