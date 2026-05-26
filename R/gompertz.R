#' Gompertz distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the Gompertz distribution.
#'
#' @details
#' The Gompertz distribution with shape \eqn{\eta > 0} and rate \eqn{b > 0} has density
#' \deqn{f(x;\,\eta,b) = b\eta\,e^{bx}\exp\!\bigl(-\eta(e^{bx}-1)\bigr), \quad x \geq 0,}
#' with CDF \eqn{F(x) = 1 - \exp(-\eta(e^{bx}-1))} and quantile function
#' \eqn{Q(p) = \log(1 - \log(1-p)/\eta)\,/\,b}.
#'
#' @param x,q vector of quantiles (non-negative)
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param eta shape parameter, must be positive
#' @param b rate parameter, must be positive
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @references
#' \url{https://en.wikipedia.org/wiki/Gompertz_distribution}
#'
#' @return
#' \code{dgompertz} gives the density, \code{pgompertz} gives the distribution function,
#' \code{qgompertz} gives the quantile function, and \code{rgompertz} generates random deviates.
#'
#' @examples
#' x <- rgompertz(1, eta = 1, b = 1)
#' d <- dgompertz(x, eta = 1, b = 1)
#' p <- pgompertz(x, eta = 1, b = 1)
#' q <- qgompertz(p, eta = 1, b = 1)
#' @name gompertz
NULL

#' @rdname gompertz
#' @export
#' @import RTMB
dgompertz <- function(x, eta = 1, b = 1, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args)
    if (any(eta <= 0)) stop("eta must be positive")
    if (any(b <= 0)) stop("b must be positive")
  }

  if (inherits(x, "simref"))
    return(dGenericSim("dgompertz", x = x, eta = eta, b = b, log = log))
  if (inherits(x, "osa"))
    return(dGenericOSA("dgompertz", x = x, eta = eta, b = b, log = log))

  logdens <- log(b) + log(eta) + b*x - eta * expm1(b*x)

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname gompertz
#' @export
#' @import RTMB
pgompertz <- function(q, eta = 1, b = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(eta <= 0)) stop("eta must be positive")
    if (any(b <= 0)) stop("b must be positive")
  }

  p <- -expm1(-eta * expm1(b*q))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}

#' @rdname gompertz
#' @export
#' @import RTMB
qgompertz <- function(p, eta = 1, b = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(eta <= 0)) stop("eta must be positive")
    if (any(b <= 0)) stop("b must be positive")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  log1p(-log1p(-p) / eta) / b
}

#' @rdname gompertz
#' @export
#' @importFrom stats rexp
rgompertz <- function(n, eta = 1, b = 1) {
  if (!ad_context()) {
    if (any(eta <= 0)) stop("eta must be positive")
    if (any(b <= 0)) stop("b must be positive")
  }

  log(1 + rexp(n, rate = eta)) / b
}
