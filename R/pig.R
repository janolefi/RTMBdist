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
#' \link[=nbinom2]{negative binomial} with \code{size = 1 / sigma}, but the
#' Poisson-inverse Gaussian has a heavier right tail and puts more mass at zero
#' for the same mean and variance. It approaches the Poisson as
#' \eqn{\sigma \to 0}.
#'
#' It is the \link[=sichel]{Sichel} distribution with \eqn{\nu = -1/2}, and is
#' computed as such; see there for how the Bessel functions are evaluated. The
#' distribution function is summed from the probability mass function, so
#' one-step-ahead (OSA) residuals via \code{RTMB::\link[RTMB]{oneStepPredict}}
#' are supported, e.g. with \code{method = "cdf"} and \code{discrete = TRUE}.
#'
#' @references
#' Dean, C., Lawless, J. F., and Willmot, G. E. (1989) A mixed Poisson-inverse-Gaussian
#' regression model. Canadian Journal of Statistics, 17, 171-181.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [sichel], [nbinom2], [invgauss]
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
dpig <- function(x, mu, sigma, log = FALSE) {

  if (!ad_context()) simulation_check(as.list(environment()))

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dpig", x = x, mu = mu, sigma = sigma, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dpig", x = x, mu = mu, sigma = sigma, log = log))
  }

  dsichel(x, mu, sigma, -0.5, log = log)
}

#' @rdname pig
#' @export
ppig <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  psichel(q, mu, sigma, -0.5, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname pig
#' @export
rpig <- function(n, mu, sigma) {
  rsichel(n, mu, sigma, -0.5)
}
