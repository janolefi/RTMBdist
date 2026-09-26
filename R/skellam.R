#' Skellam distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the Skellam distribution.
#'
#' The Skellam distribution is the distribution of the difference of two Poisson random variables. Specifically, if \eqn{X_1 \sim \text{Pois}(\mu_1)} and \eqn{X_2 \sim \text{Pois}(\mu_2)}, then \eqn{X_1 - X_2 \sim \text{Skellam}(\mu_1, \mu_2)}.
#'
#' @details
#' This implementation of \code{dskellam} allows for automatic differentiation with \code{RTMB}.
#'
#' The distribution function has no closed form, and as the support is unbounded in both
#' directions, summing the probability mass function would need a range that depends on
#' the parameters. Instead, \code{pskellam} uses the identity
#' \deqn{P(X \le q;\, \mu_1, \mu_2) = \int_{\mu_1}^\infty P(X = q;\, t, \mu_2)\, dt,}
#' which follows from the relation between the Poisson and gamma distribution functions,
#' and integrates it numerically with the AD-compatible
#' \code{\link[RTMB:ADintegrate]{integrate}} of \code{RTMB}. By symmetry, the upper tail
#' is \eqn{P(X > q;\, \mu_1, \mu_2) = P(X \le -q - 1;\, \mu_2, \mu_1)}, and whichever
#' tail is smaller is integrated directly, which keeps small tail probabilities accurate.
#' \code{pskellam} is AD-compatible in \code{mu1} and \code{mu2}, while \code{q} must be
#' numeric data. This is also what one-step-ahead (OSA) residuals via
#' \code{RTMB::\link[RTMB]{oneStepPredict}} need, so these are supported, e.g. with
#' \code{method = "cdf"} and \code{discrete = TRUE}.
#'
#' @param x integer vector of counts
#' @param q vector of quantiles
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le q]}, otherwise \eqn{P[X > q]}.
#' @param log.p logical; if \code{TRUE}, probabilities are returned on the log scale.
#' @param n number of random values to return.
#' @param mu1,mu2 Poisson means
#' @param log logical; return log-density if TRUE
#'
#' @return
#' \code{dskellam} gives the probability mass function, \code{pskellam} gives the distribution function, and \code{rskellam} generates random deviates.
#'
#' @examples
#' x <- rskellam(1, 2, 3)
#' d <- dskellam(x, 2, 3)
#' p <- pskellam(x, 2, 3)
#' @name skellam
NULL
#' @rdname skellam
#' @export
#' @import RTMB
dskellam <- function(x, mu1, mu2, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure lambda, phi > 0
    if (any(mu1 <= 0)) stop("mu1 must be > 0")
    if (any(mu2 <= 0)) stop("mu2 must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dskellam", x = x, mu1 = mu1, mu2 = mu2, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dskellam", x = x, mu1 = mu1, mu2 = mu2, log = log))
  }

  val <- 2 * sqrt(mu1 * mu2)
  logI <- log(besselI(val, x, expon.scaled = TRUE)) + val

  logdens <- -mu1 - mu2 + (x / 2) * (log(mu1) - log(mu2)) + logI

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname skellam
#' @export
pskellam <- function(q, mu1, mu2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(mu1 <= 0)) stop("mu1 must be > 0")
    if (any(mu2 <= 0)) stop("mu2 must be > 0")
  }
  if (inherits(q, "advector")) stop("q must be numeric data, not an AD variable.")

  n <- max(length(q), length(mu1), length(mu2))
  rec <- function(x) x[(seq_len(n) - 1) %% length(x) + 1]
  q <- floor(rec(q))
  mu1 <- rec(mu1)
  mu2 <- rec(mu2)

  # The lower tail is P(X <= q; mu1, mu2) = int_{mu1}^Inf dskellam(q, t, mu2) dt, and by
  # symmetry the upper tail is P(X > q; mu1, mu2) = P(X <= -q - 1; mu2, mu1). The smaller one
  # is integrated directly: the lower tail below the mean mu1 - mu2, the upper one above.
  # Only the values of the parameters at tape time decide this, so the choice stays fixed
  # when the tape is re-evaluated, which leaves the result correct and only affects accuracy.
  lower <- q < value_of(mu1) - value_of(mu2)
  h <- sqrt(value_of(mu1) + value_of(mu2)) # unit of the integration, of the order of the sd

  # Log of int_{m1}^Inf dskellam(k, t, m2) dt, integrated as h int_0^Inf over t = m1 + h v,
  # so that the limits are fixed: RTMB's AD integrate() can give a zero derivative with
  # respect to a finite limit, and m1 is a parameter. The integrand is taken relative to its
  # value at v = 0, where it is largest, so that the integral is of order one far out in
  # the tails as well, which avoids underflow (subnormal values give NaN derivatives) and
  # keeps log.p accurate. Each element gets its own integrate() call, as the second
  # derivatives through RTMB's Vectorize() are wrong for integrands like this (RTMB 2.0).
  log_tail <- function(k, m1, m2, h) {
    lp0 <- dskellam(k, m1, m2, log = TRUE)
    r <- function(v) exp(dskellam(k, m1 + h * v, m2, log = TRUE) - lp0)
    lp0 + log(h * integrate(r, 0, Inf, rel.tol = 1e-10, abs.tol = 0)$value)
  }
  # the directly integrated tail, from P(X = k; m1, m2) with k = q or -q - 1
  # (swapped by 0/1 weights rather than by assignment, which would fail when only one of
  # mu1 and mu2 is an AD variable)
  k <- ifelse(lower, q, -q - 1)
  m1 <- lower * mu1 + (1 - lower) * mu2
  m2 <- lower * mu2 + (1 - lower) * mu1

  # This tail is 0 for infinite q. It is also treated as 0 where the probability mass
  # function itself underflows at tape time, i.e. is below about 1e-308, as the ratio in
  # the integrand cannot be formed there.
  ls <- rep(-Inf, n)
  if (ad_context()) ls <- advector(ls)
  ok <- which(is.finite(q))
  ok <- ok[is.finite(value_of(dskellam(k[ok], m1[ok], m2[ok], log = TRUE)))]
  if (length(ok)) ls[ok] <- do.call(c, lapply(ok, function(i) log_tail(k[i], m1[i], m2[i], h[i])))

  flip <- lower != lower.tail
  p <- if (log.p) ls else exp(ls)
  if (any(flip, na.rm = TRUE)) {
    i <- which(flip)
    p[i] <- if (log.p) log1p(-exp(ls[i])) else 1 - exp(ls[i])
  }
  if (anyNA(q)) p[is.na(q)] <- NA
  p
}

#' @rdname skellam
#' @export
rskellam <- function(n, mu1, mu2) {

  if (any(mu1 <= 0)) stop("mu1 must be > 0")
  if (any(mu2 <= 0)) stop("mu2 must be > 0")

  x1 <- rpois(n, mu1)
  x2 <- rpois(n, mu2)

  return(x1 - x2)
}

