#' Wrapped normal distribution
#'
#' Density, distribution function, and random generation for the wrapped normal
#' distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' The wrapped normal is the distribution of \eqn{Y \bmod 2\pi} for a normal
#' \eqn{Y} with mean \eqn{\mu} and standard deviation \eqn{\sigma}, with density
#' \deqn{f(x;\,\mu,\sigma) = \sum_{k=-\infty}^{\infty} \frac{1}{\sigma}\,\phi\!\left(\frac{x - \mu + 2\pi k}{\sigma}\right)
#'       = \frac{1}{2\pi}\left(1 + 2\sum_{k=1}^{\infty} e^{-k^2\sigma^2/2} \cos(k(x - \mu))\right).}
#' For \eqn{\sigma < 2} the first sum is used with \eqn{|k| \le 3}, otherwise the
#' second with \eqn{k \le 4}; the terms left out are below \eqn{10^{-20}}.
#'
#' The distribution function cuts the circle open at \eqn{\mu - \pi} by default,
#' or at \code{from}, exactly as \code{\link{pwrpcauchy}} does, where the choice
#' of origin and the one-step-ahead (OSA) residuals are explained. As there, OSA
#' residuals via \code{RTMB::\link[RTMB]{oneStepPredict}} with
#' \code{method = "cdf"} cut the circle at the fixed origin \eqn{-\pi}.
#'
#' @seealso [wrpcauchy], [vm], [card]; [cjw()] for circular-linear copulas.
#'
#' @param x,q vector of angles in radians.
#' @param mu mean direction in radians.
#' @param sigma standard deviation of the normal distribution before wrapping, must be positive.
#' @param from origin, in radians, at which the circle is cut open for the distribution function. If \code{NULL} (default), it is set to \code{mu - pi}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param n number of random values to return.
#'
#' @return \code{dwrpnorm} gives the density, \code{pwrpnorm} gives the distribution function, and \code{rwrpnorm} generates random deviates in \eqn{[-\pi, \pi)}.
#'
#' @examples
#' set.seed(1)
#' x <- rwrpnorm(10, 0, 0.8)
#' d <- dwrpnorm(x, 0, 0.8)
#' p <- pwrpnorm(x, 0, 0.8)
#' @name wrpnorm
NULL

#' @rdname wrpnorm
#' @export
dwrpnorm <- function(x, mu = 0, sigma, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(sigma <= 0)) stop("sigma must be positive.")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dwrpnorm", x = x, mu = mu, sigma = sigma, log = log))
  }
  if (inherits(x, "osa")) {
    # the circle is cut at the fixed origin -pi, as for dwrpcauchy
    return(dGenericOSA("dwrpnorm_osa", x = x, mu = mu, sigma = sigma, log = log))
  }

  u <- wrap_angle(x - mu)
  small <- smaller(sigma, 2)
  s1 <- sigma - (1 - small) * (sigma - 2) # sigma, at most 2
  s2 <- sigma + small * (2 - sigma) # sigma, at least 2

  # small sigma: the windings of the normal density closest to u
  ld1 <- dnorm(u, 0, s1, log = TRUE)
  for (k in c(-3:-1, 1:3)) ld1 <- logspace_add(ld1, dnorm(u + 2 * pi * k, 0, s1, log = TRUE))

  # large sigma: Fourier series
  f2 <- 1
  for (k in 1:4) f2 <- f2 + 2 * exp(-k * k * s2 * s2 / 2) * cos(k * u)
  ld2 <- log(f2 / (2 * pi))

  logdens <- small * ld1 + (1 - small) * ld2

  if (log) return(logdens)
  return(exp(logdens))
}

# density and distribution function behind the OSA residuals of dwrpnorm, with the
# circle cut at -pi, found by RTMB's dGenericOSA via their names
dwrpnorm_osa <- function(x, mu, sigma, log = FALSE) dwrpnorm(x, mu, sigma, log = log)
pwrpnorm_osa <- function(q, mu, sigma) pwrpnorm(q, mu, sigma, from = -pi)

#' @rdname wrpnorm
#' @export
pwrpnorm <- function(q, mu = 0, sigma, from = NULL, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be positive.")
  }

  small <- smaller(sigma, 2)
  s1 <- sigma - (1 - small) * (sigma - 2)
  s2 <- sigma + small * (2 - sigma)

  # probability of [mu - pi, mu + u), from the same two series as the density
  F0 <- function(u) {
    F1 <- 0
    for (k in -3:3) F1 <- F1 + pnorm((u + 2 * pi * k) / s1) - pnorm((2 * k - 1) * pi / s1)
    F2 <- (u + pi) / (2 * pi)
    for (k in 1:4) F2 <- F2 + exp(-k * k * s2 * s2 / 2) * sin(k * u) / (pi * k)
    small * F1 + (1 - small) * F2
  }

  pcircular(F0, q, mu, from, lower.tail, log.p)
}

#' @rdname wrpnorm
#' @export
#' @importFrom stats rnorm
rwrpnorm <- function(n, mu = 0, sigma) {

  if (any(sigma <= 0)) stop("sigma must be positive.")

  wrap_angle(rnorm(n, mu, sigma))
}
