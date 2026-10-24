#' Poisson-inverse Gaussian distribution
#'
#' Probability mass function, distribution function and random generation for
#' the Poisson-inverse Gaussian distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' The Poisson-inverse Gaussian is a mixed Poisson distribution whose mean is
#' drawn from an inverse Gaussian distribution with mean \eqn{\mu} and variance
#' \eqn{\sigma \mu^2}. Its probability mass function is
#' \deqn{P(X = x) = \left(\frac{2\alpha}{\pi}\right)^{1/2}
#'       \frac{\mu^x e^{1/\sigma} K_{x - 1/2}(\alpha)}{(\alpha\sigma)^x\, x!},
#'       \qquad \alpha^2 = \frac{1}{\sigma^2} + \frac{2\mu}{\sigma},
#'       \quad x = 0, 1, 2, \ldots}
#' where \eqn{K_\nu} is the modified Bessel function of the second kind. The
#' parameterisation follows the \code{PIG} family of the \code{gamlss.dist}
#' package, with
#' \deqn{E(X) = \mu, \qquad \mathrm{Var}(X) = \mu + \sigma \mu^2.}
#' This is the same mean-variance relationship as the
#' \link[=nbinom2]{negative binomial} with \code{phi = 1 / sigma}, but the
#' Poisson-inverse Gaussian has a heavier right tail and puts more mass at zero
#' for the same mean and variance. It approaches the Poisson as
#' \eqn{\sigma \to 0}.
#'
#' The order of the Bessel function depends only on the data, so no derivative
#' with respect to it is needed. For counts below 30 it is evaluated with
#' \code{RTMB}'s \code{besselK}. Above that, \code{besselK} overflows when
#' \eqn{\alpha} is small, so its uniform asymptotic expansion for large orders
#' is used instead, which is accurate to about \eqn{10^{-11}} on the log scale
#' there.
#'
#' The distribution function has no closed form and is computed by summing the
#' probability mass function over \eqn{0, \ldots, q}, so its cost grows with the
#' largest \code{q}. It is AD-compatible in the parameters, while \code{q} must be
#' numeric data. This is also what one-step-ahead (OSA) residuals via
#' \code{RTMB::\link[RTMB]{oneStepPredict}} need, so these are supported, e.g. with
#' \code{method = "cdf"} and \code{discrete = TRUE}.
#'
#' @references
#' Dean, C., Lawless, J. F., and Willmot, G. E. (1989) A mixed Poisson-inverse-Gaussian
#' regression model. Canadian Journal of Statistics, 17, 171-181.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' Olver, F. W. J. et al. (eds.) NIST Digital Library of Mathematical Functions,
#' equation 10.41.4, https://dlmf.nist.gov/10.41.
#'
#' @seealso [nbinom2], [invgauss], [bnbinom2]
#'
#' @param x vector of non-negative counts.
#' @param q vector of quantiles.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le q]}, otherwise \eqn{P[X > q]}.
#' @param log.p logical; if \code{TRUE}, probabilities are returned on the log scale.
#' @param n number of random values to return (for \code{rpig}).
#' @param mu mean parameter, must be positive.
#' @param sigma dispersion parameter, must be positive.
#' @param log logical; if \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return
#' \code{dpig} gives the probability mass function, \code{ppig} gives the distribution function, and \code{rpig} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rpig(5, mu = 4, sigma = 0.5)
#' d <- dpig(x, mu = 4, sigma = 0.5)
#' p <- ppig(x, mu = 4, sigma = 0.5)
#' @name pig
NULL

#' @rdname pig
#' @export
#' @import RTMB
dpig <- function(x, mu, sigma, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(mu <= 0)) stop("mu must be positive.")
    if (any(sigma <= 0)) stop("sigma must be positive.")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dpig", x = x, mu = mu, sigma = sigma, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dpig", x = x, mu = mu, sigma = sigma, log = log))
  }

  # clamped below zero so that lgamma stays away from its poles; the density
  # itself is set to zero there by the indicator
  xc <- pmax.ad(x, 0)

  # alpha * sigma = sqrt(1 + 2 mu sigma). The factor exp(1 / sigma) of the mass
  # function is combined with the exp(-alpha) taken out of the scaled Bessel
  # function, and 1 / sigma - alpha is written without the cancellation
  s <- sqrt(1 + 2 * mu * sigma)
  alpha <- s / sigma

  logdens <- 0.5 * log(2 * alpha / pi) + xc * (log(mu) - log(s)) - lgamma(xc + 1) -
    2 * mu / (1 + s) + log_besselK_scaled(alpha, xc - 0.5) +
    log(1 - smaller(x, 0)) # zero below the support

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname pig
#' @export
ppig <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(mu <= 0)) stop("mu must be positive.")
    if (any(sigma <= 0)) stop("sigma must be positive.")
  }

  discrete_cdf(dpig, q, list(mu = mu, sigma = sigma),
               lower.tail = lower.tail, log.p = log.p)
}

#' @rdname pig
#' @export
#' @importFrom stats rpois
rpig <- function(n, mu, sigma) {

  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 0)) stop("sigma must be positive.")

  n <- ceiling(n)

  # draw the Poisson mean from an inverse Gaussian with mean mu and variance sigma * mu^2,
  # as mu times one with mean 1 and shape 1 / sigma
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  lambda <- mu * statmod::rinvgauss(n, mean = 1, shape = 1 / sigma)

  rpois(n, lambda)
}

# Log of the exponentially scaled modified Bessel function of the second kind,
# log(besselK(x, nu, expon.scaled = TRUE)), for numeric orders nu, AD-compatible in x.
# besselK() overflows for large orders and small x, so from |nu| >= 29.5 on the
# uniform asymptotic (Debye) expansion for large orders, DLMF 10.41.4, is used with
# its first six correction terms, which keeps the absolute error on the log scale
# below about 1e-11 there. Both branches are evaluated with the order clamped to their
# range, so that neither produces Inf, and combined with indicators.
log_besselK_scaled <- function(x, nu, cut = 29.5) {
  nu <- abs(nu) # K_{-nu} = K_nu
  big <- 1 - smaller(nu, cut)
  nu_s <- nu - big * (nu - cut) # clamped to <= cut
  nu_b <- nu + (1 - big) * (cut - nu) # clamped to >= cut

  # small orders
  small_part <- log(besselK(x, nu_s, expon.scaled = TRUE))

  # large orders: K_nu(nu z) ~ sqrt(pi / (2 nu)) exp(-nu eta) / (1 + z^2)^(1/4) *
  # sum_k (-1)^k u_k(t) / nu^k, with t = 1 / sqrt(1 + z^2)
  z <- x / nu_b
  r <- sqrt(1 + z * z)
  t <- 1 / r
  t2 <- t * t
  eta <- r + log(z / (1 + r))
  poly <- function(coef) { # polynomial in t^2 by Horner's scheme
    res <- 0
    for (a in rev(coef)) res <- res * t2 + a
    res
  }
  u1 <- t * poly(c(3, -5)) / 24
  u2 <- t2 * poly(c(81, -462, 385)) / 1152
  u3 <- t * t2 * poly(c(30375, -369603, 765765, -425425)) / 414720
  u4 <- t2 * t2 * poly(c(4465125, -94121676, 349922430, -446185740, 185910725)) / 39813120
  u5 <- t * t2 * t2 * poly(c(1519035525, -49286948607, 284499769554, -614135872350,
                             566098157625, -188699385875)) / 6688604160
  u6 <- t2 * t2 * t2 * poly(c(2757049477875, -127577298354750, 1050760774457901,
                              -3369032068261860, 5104696716244125, -3685299006138750,
                              1023694168371875)) / 4815794995200
  v <- 1 / nu_b
  series <- 1 + v * (-u1 + v * (u2 + v * (-u3 + v * (u4 + v * (-u5 + v * u6)))))
  big_part <- 0.5 * log(pi / (2 * nu_b)) - nu_b * eta - 0.5 * log(r) + log(series) + x

  big * big_part + (1 - big) * small_part
}
