#' Gumbel distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the Gumbel distribution.
#'
#' @details
#' This implementation of \code{dgumbel} allows for automatic differentiation with \code{RTMB}.
#'
#' \deqn{f(x;\,\mu,\sigma) = \frac{1}{\sigma}\exp\!\bigl(-(z + e^{-z})\bigr), \quad z = \frac{x-\mu}{\sigma}.}
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param location location parameter
#' @param scale scale parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dgumbel} gives the density, \code{pgumbel} gives the distribution function, \code{qgumbel} gives the quantile function, and \code{rgumbel} generates random deviates.
#'
#' @examples
#' x <- rgumbel(1, 0.5, 2)
#' d <- dgumbel(x, 0.5, 2)
#' p <- pgumbel(x, 0.5, 2)
#' q <- qgumbel(p, 0.5, 2)
#' @name gumbel
NULL

#' @rdname gumbel
#' @export
dgumbel <- function(x, location = 0, scale = 1, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if(scale <= 0) stop("Scale parameter must be positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dgumbel", x=x, location=location, scale=scale, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dgumbel", x=x, location=location, scale=scale, log=log))
  }

  z <- (x - location) / scale
  logdens <- - z - exp(-z) - log(scale)

  if(log) return(logdens)

  return(exp(logdens))
}
#' @rdname gumbel
#' @export
pgumbel <- function(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if(!ad_context()) {
    if(scale <= 0) stop("Scale parameter must be positive.")
  }

  z <- (q - location) / scale

  p <- exp(-exp(-z))

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname gumbel
#' @export
qgumbel <- function(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if(!ad_context()) {
    if(scale <= 0) stop("Scale parameter must be positive.")
  }

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p

  location - scale * log(-log(p))
}
#' @rdname gumbel
#' @export
#' @importFrom stats runif
rgumbel <- function(n, location = 0, scale = 1) {
  if(!ad_context()) {
    if(scale <= 0) stop("Scale parameter must be positive.")
  }

  u <- runif(n)
  qgumbel(u, location = location, scale = scale)
}


