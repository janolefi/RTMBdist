#' Kumaraswamy distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the Kumaraswamy distribution.
#'
#' @param x,q vector of quantiles in \eqn{(0,1)}
#' @param p vector of probabilities
#' @param a,b positive shape parameters
#' @param n number of random values to return.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @returns
#' \code{dkumar} gives the density, \code{pkumar} gives the distribution function, \code{qkumar} gives the quantile function, and \code{rkumar} generates random deviates.
#'
#' @examples
#' x <- rkumar(1, a = 1, b = 2)
#' d <- dkumar(x, a = 1, b = 2)
#' p <- pkumar(x, a = 1, b = 2)
#' q <- qkumar(p, a = 1, b = 2)
#' @name kumar
NULL

#' @rdname kumar
#' @export
#' @import RTMB
dkumar <- function(x, a, b, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(a <= 0))  stop("a must be positive")
    if (any(b <= 0))  stop("b must be positive")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dkumar", x=x, a=a, b=b, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dkumar", x=x, a=a, b=b, log=log))
  }

  logdens <- log(a) + log(b) + (a-1) * log(x) + (b-1) * log1p(-x^a)

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname kumar
#' @export
#' @import RTMB
pkumar <- function(q, a, b, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(a <= 0))  stop("a must be positive")
    if (any(b <= 0))  stop("b must be positive")
  }

  p <- -expm1(b * log1p(-q^a)) # more stable, avoid double power

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}

#' @rdname kumar
#' @export
#' @import RTMB
qkumar <- function(p, a, b, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(a <= 0))  stop("a must be positive")
    if (any(b <= 0))  stop("b must be positive")
  }

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p

  exp(log1p(- (1-p)^(1/b)) / a) # more stable, avoid double power
}

#' @rdname kumar
#' @export
#' @importFrom stats runif
rkumar <- function(n, a, b) {

  args <- as.list(environment())
  simulation_check(args) # informative error message if likelihood in wrong order
  if (any(a <= 0))  stop("a must be positive")
  if (any(b <= 0))  stop("b must be positive")

  n <- ceiling(n)
  p <- runif(n)
  qkumar(p, a=a, b=b)
}
