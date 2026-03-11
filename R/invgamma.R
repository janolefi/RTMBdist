#' Inverse Gamma distribution
#'
#' Density, distribution function, and random generation for
#' the inverse Gamma distribution.
#'
#' @details
#' This implementation of \code{dinvgamma}, \code{pinvgamma}, and \code{qinvgamma} allows for automatic differentiation with \code{RTMB}.
#'
#' If \eqn{X \sim \Gamma(\alpha, \beta)}, then \eqn{1/X \sim \text{InvGamma}(\alpha, \beta)}.
#'
#' @param x,q vector of quantiles, must be positive.
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param shape,rate,scale positive parameters of corresponding gamma distribution
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dinvgamma} gives the density, \code{pinvgamma} gives the distribution function, \code{qinvgamma} gives the quantile function, and \code{rinvgamma} generates random deviates.
#'
#' @examples
#' x <- rinvgamma(1, 1, 0.5)
#' d <- dinvgamma(x, 1, 0.5)
#' p <- pinvgamma(x, 1, 0.5)
#' q <- qinvgamma(p, 1, 0.5)
#' @name invgamma
NULL

#' @rdname invgamma
#' @export
#' @import RTMB
dinvgamma <- function(x, shape, rate, scale = 1/rate, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure shape, rate, scale > 0
    if (any(shape <= 0)) stop("shape must be strictly positive.")
    if (any(rate <= 0)) stop("rate must be strictly positive.")
    if (any(scale <= 0)) stop("scale must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dinvgamma", x=x, shape=shape, rate=rate, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dinvgamma", x=x, shape=shape, rate=rate, log=log))
  }

  logC <- shape * log(scale) - lgamma(shape)
  logdens <- - (shape + 1) * log(x) - scale / x + logC

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname invgamma
#' @export
#' @import RTMB
pinvgamma <- function(q, shape, rate, scale = 1/rate, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure shape, rate, scale > 0
    if (any(shape <= 0)) stop("shape must be strictly positive.")
    if (any(rate <= 0)) stop("rate must be strictly positive.")
    if (any(scale <= 0)) stop("scale must be strictly positive.")
  }

  p <- reggamma(shape, scale / q)

  if (!lower.tail) p <- 1 - p
  if (log.p) return(log(p))
  return(p)
}

#' @rdname invgamma
#' @export
#' @import RTMB
qinvgamma <- function(p, shape, rate, scale = 1/rate, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure shape, rate, scale > 0
    if (any(shape <= 0)) stop("shape must be strictly positive.")
    if (any(rate <= 0)) stop("rate must be strictly positive.")
    if (any(scale <= 0)) stop("scale must be strictly positive.")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  q <- RTMB::qgamma(1 - p, shape = shape, scale = 1)
  scale / q
}

#' @rdname invgamma
#' @export
#' @importFrom stats rgamma
rinvgamma <- function(n, shape, rate, scale = 1/rate) {

  # ensure shape, rate, scale > 0
  if (any(shape <= 0)) stop("shape must be strictly positive.")
  if (any(rate <= 0)) stop("rate must be strictly positive.")
  if (any(scale <= 0)) stop("scale must be strictly positive.")

  x <- stats::rgamma(n, shape=shape, rate=rate)
  1 / x
}
