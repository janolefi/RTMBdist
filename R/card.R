#' Cardioid distribution
#'
#' Density, distribution function, and random generation for the cardioid
#' distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' The cardioid distribution has density
#' \deqn{f(x;\,\mu,\rho) = \frac{1 + 2\rho\cos(x - \mu)}{2\pi}, \qquad 0 \le \rho \le \frac{1}{2},}
#' and is the circular uniform for \eqn{\rho = 0}. With the circle cut open at
#' \eqn{\mu - \pi}, its distribution function is
#' \deqn{F(q) = \frac{q - \mu + \pi}{2\pi} + \frac{\rho}{\pi}\sin(q - \mu), \quad \mu - \pi \le q < \mu + \pi.}
#' The origin can be changed with \code{from}, exactly as for
#' \code{\link{pwrpcauchy}}, where the choice of origin and the one-step-ahead
#' (OSA) residuals are explained. As there, OSA residuals via
#' \code{RTMB::\link[RTMB]{oneStepPredict}} with \code{method = "cdf"} cut the
#' circle at the fixed origin \eqn{-\pi}.
#'
#' @references
#' Jammalamadaka, S. R. and SenGupta, A. (2001) Topics in Circular Statistics. World Scientific.
#'
#' @seealso [wrpcauchy], [wrpnorm], [vm]
#'
#' @param x,q vector of angles in radians.
#' @param mu mean direction in radians.
#' @param rho concentration parameter, must be in the interval from 0 to 1/2.
#' @param from origin, in radians, at which the circle is cut open for the distribution function. If \code{NULL} (default), it is set to \code{mu - pi}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param n number of random values to return.
#'
#' @return \code{dcard} gives the density, \code{pcard} gives the distribution function, and \code{rcard} generates random deviates in \eqn{[-\pi, \pi)}.
#'
#' @examples
#' set.seed(1)
#' x <- rcard(10, 0, 0.3)
#' d <- dcard(x, 0, 0.3)
#' p <- pcard(x, 0, 0.3)
#' @name card
NULL

#' @rdname card
#' @export
dcard <- function(x, mu = 0, rho, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(rho < 0) || any(rho > 0.5)) stop("rho must be in the interval [0, 0.5].")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dcard", x = x, mu = mu, rho = rho, log = log))
  }
  if (inherits(x, "osa")) {
    # the circle is cut at the fixed origin -pi, as for dwrpcauchy
    return(dGenericOSA("dcard_osa", x = x, mu = mu, rho = rho, log = log))
  }

  logdens <- log1p(2 * rho * cos(x - mu)) - log(2 * pi)

  if (log) return(logdens)
  return(exp(logdens))
}

# density and distribution function behind the OSA residuals of dcard, with the
# circle cut at -pi, found by RTMB's dGenericOSA via their names
dcard_osa <- function(x, mu, rho, log = FALSE) dcard(x, mu, rho, log = log)
pcard_osa <- function(q, mu, rho) pcard(q, mu, rho, from = -pi)

#' @rdname card
#' @export
pcard <- function(q, mu = 0, rho, from = NULL, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(rho < 0) || any(rho > 0.5)) stop("rho must be in the interval [0, 0.5].")
  }

  F0 <- function(u) (u + pi) / (2 * pi) + rho * sin(u) / pi

  pcircular(F0, q, mu, from, lower.tail, log.p)
}

#' @rdname card
#' @export
#' @importFrom stats runif
rcard <- function(n, mu = 0, rho) {

  if (any(rho < 0) || any(rho > 0.5)) stop("rho must be in the interval [0, 0.5].")

  n <- ceiling(n)
  mu <- rep_len(mu, n)
  rho <- rep_len(rho, n)

  # rejection sampling of the angle u from mu, with a uniform proposal
  x <- numeric(n)
  todo <- seq_len(n)
  while (length(todo)) {
    u <- runif(length(todo), -pi, pi)
    keep <- runif(length(todo)) * (1 + 2 * rho[todo]) <= 1 + 2 * rho[todo] * cos(u)
    x[todo[keep]] <- wrap_angle(mu[todo[keep]] + u[keep])
    todo <- todo[!keep]
  }
  x
}
