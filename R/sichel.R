#' Sichel distribution
#'
#' Probability mass function, distribution function and random generation for
#' the Sichel distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' The Sichel distribution is a mixed Poisson distribution whose mean follows a
#' generalised inverse Gaussian distribution. In the parameterisation of the
#' \code{SICHEL} family of the \code{gamlss.dist} package, \eqn{\mu} is the mean
#' and the probability mass function is
#' \deqn{P(X = x) = \frac{(\mu/c)^x\, K_{x + \nu}(\alpha)}{x!\,(\alpha\sigma)^{x + \nu} K_\nu(1/\sigma)},
#'       \qquad \alpha^2 = \frac{1}{\sigma^2} + \frac{2\mu}{c\,\sigma},
#'       \quad c = \frac{K_{\nu + 1}(1/\sigma)}{K_\nu(1/\sigma)},}
#' for \eqn{x = 0, 1, 2, \ldots}, where \eqn{K_\nu} is the modified Bessel
#' function of the second kind. The variance is
#' \deqn{\mathrm{Var}(X) = \mu + \mu^2\left(\frac{2\sigma(\nu + 1)}{c} + \frac{1}{c^2} - 1\right).}
#' For \eqn{\nu = -1/2} this is the \link[=pig]{Poisson-inverse Gaussian}
#' distribution, and as \eqn{\sigma \to \infty} with \eqn{\nu > 0} it approaches
#' the \link[=nbinom2]{negative binomial}.
#'
#' The Bessel functions are evaluated with \code{RTMB}'s \code{besselK} for orders
#' up to 50. For larger orders, \code{besselK} overflows when its argument is
#' small, so its uniform asymptotic expansion for large orders is used instead,
#' which is accurate to about \eqn{10^{-10}} on the log scale there. Both are
#' differentiable in the order, so \eqn{\nu} can be estimated.
#'
#' The distribution function has no closed form and is computed by summing the
#' probability mass function over \eqn{0, \ldots, q}, so its cost grows with the
#' largest \code{q}. It is AD-compatible in the parameters, while \code{q} must be
#' numeric data. This is also what one-step-ahead (OSA) residuals via
#' \code{RTMB::\link[RTMB]{oneStepPredict}} need, so these are supported, e.g. with
#' \code{method = "cdf"} and \code{discrete = TRUE}.
#'
#' @references
#' Sichel, H. S. (1982) Repeat-buying and the generalized inverse Gaussian-Poisson
#' distribution. Journal of the Royal Statistical Society, Series C, 31, 193-204.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' Olver, F. W. J. et al. (eds.) NIST Digital Library of Mathematical Functions,
#' equation 10.41.4, https://dlmf.nist.gov/10.41.
#'
#' @seealso [pig], [nbinom2], [bnbinom2]
#'
#' @param x vector of non-negative counts.
#' @param q vector of quantiles.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le q]}, otherwise \eqn{P[X > q]}.
#' @param log.p logical; if \code{TRUE}, probabilities are returned on the log scale.
#' @param n number of random values to return (for \code{rsichel}).
#' @param mu mean parameter, must be positive.
#' @param sigma dispersion parameter, must be positive.
#' @param nu shape parameter, any real number.
#' @param log logical; if \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return
#' \code{dsichel} gives the probability mass function, \code{psichel} gives the distribution function, and \code{rsichel} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rsichel(5, mu = 4, sigma = 0.5, nu = 1)
#' d <- dsichel(x, mu = 4, sigma = 0.5, nu = 1)
#' p <- psichel(x, mu = 4, sigma = 0.5, nu = 1)
#' @name sichel
NULL

#' @rdname sichel
#' @export
#' @import RTMB
dsichel <- function(x, mu, sigma, nu, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    check_sichel(mu, sigma)
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dsichel", x = x, mu = mu, sigma = sigma, nu = nu, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dsichel", x = x, mu = mu, sigma = sigma, nu = nu, log = log))
  }

  # clamped below zero so that lgamma stays away from its poles; the mass function
  # itself is set to zero there by the indicator
  xc <- pmax.ad(x, 0)

  # The Bessel functions are exponentially scaled. With s = alpha * sigma, the factors
  # exp(1 / sigma - alpha) they leave combine to exp(-2 mu / (c (1 + s))), which avoids
  # the cancellation in 1 / sigma - alpha
  lK <- log_besselK_scaled(1 / sigma, nu)
  c <- exp(log_besselK_scaled(1 / sigma, nu + 1) - lK)
  s <- sqrt(1 + 2 * mu * sigma / c)

  logdens <- xc * log(mu / c) - (xc + nu) * log(s) - lgamma(xc + 1) -
    2 * mu / (c * (1 + s)) + log_besselK_scaled(s / sigma, xc + nu) - lK +
    log(1 - smaller(x, 0)) # zero below the support

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname sichel
#' @export
psichel <- function(q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) check_sichel(mu, sigma)

  discrete_cdf(dsichel, q, list(mu = mu, sigma = sigma, nu = nu),
               lower.tail = lower.tail, log.p = log.p)
}

#' @rdname sichel
#' @export
#' @importFrom stats runif
rsichel <- function(n, mu, sigma, nu) {

  check_sichel(mu, sigma)

  n <- ceiling(n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  nu <- rep_len(nu, n)

  # inversion: the smallest x with P(X <= x) >= u, found by summing the mass function.
  # A draw also stops once its terms underflow, so rounding cannot make the loop endless
  u <- runif(n)
  x <- numeric(n)
  F <- dsichel(x, mu, sigma, nu)
  todo <- which(F < u)
  while (length(todo)) {
    x[todo] <- x[todo] + 1
    d <- dsichel(x[todo], mu[todo], sigma[todo], nu[todo])
    F[todo] <- F[todo] + d
    todo <- todo[F[todo] < u[todo] & d > 0]
  }
  x
}

check_sichel <- function(mu, sigma) {
  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 0)) stop("sigma must be positive.")
}

# log(besselK(x, nu, expon.scaled = TRUE)), AD-compatible in x and nu.
# besselK() overflows for large orders and small x, so from |nu| = 50 on the uniform
# asymptotic expansion for large orders (DLMF 10.41.4) is used, with four correction
# terms; its error on the log scale is below 1e-10 there. Each branch gets the order
# clamped to its own range, so that neither produces Inf.
log_besselK_scaled <- function(x, nu, cut = 50) {
  nu <- abs(nu) # K_{-nu} = K_nu
  big <- 1 - smaller(nu, cut)
  nu_s <- nu - big * (nu - cut) # at most cut
  nu_b <- nu + (1 - big) * (cut - nu) # at least cut

  small_part <- log(besselK(x, nu_s, expon.scaled = TRUE))

  z <- x / nu_b
  r <- sqrt(1 + z * z)
  t <- 1 / r
  u1 <- (3 * t - 5 * t^3) / 24
  u2 <- (81 * t^2 - 462 * t^4 + 385 * t^6) / 1152
  u3 <- (30375 * t^3 - 369603 * t^5 + 765765 * t^7 - 425425 * t^9) / 414720
  u4 <- (4465125 * t^4 - 94121676 * t^6 + 349922430 * t^8 - 446185740 * t^10 +
           185910725 * t^12) / 39813120
  series <- 1 - u1 / nu_b + u2 / nu_b^2 - u3 / nu_b^3 + u4 / nu_b^4
  big_part <- 0.5 * log(pi / (2 * nu_b)) - nu_b * (r + log(z / (1 + r))) -
    0.5 * log(r) + log(series) + x

  big * big_part + (1 - big) * small_part
}
