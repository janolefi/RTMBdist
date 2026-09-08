#' Reparameterised beta-negative binomial distribution
#'
#' Probability mass function and random generation for the beta-negative
#' binomial distribution reparameterised in terms of its mean.
#'
#' @details
#' \code{dbnbinom2} allows for automatic differentiation with \code{RTMB}.
#' The parameterisation follows the \code{BNB} family of the \code{gamlss.dist}
#' package, in which \eqn{\mu} is exactly the mean.
#'
#' Writing \eqn{r}, \eqn{a} and \eqn{b} for the arguments \code{size},
#' \code{shape1} and \code{shape2} of \code{\link{dbnbinom}}, the
#' reparameterisation is
#' \deqn{a = \frac{1}{\sigma} + 1, \qquad b = \frac{\mu\nu}{\sigma}, \qquad r = \frac{1}{\nu}.}
#'
#' This gives \eqn{E(X) = \mu} for every admissible \eqn{\sigma} and \eqn{\nu},
#' where the original parameterisation needs \eqn{a > 1} for a mean to exist at
#' all. The variance is
#' \deqn{\mathrm{Var}(X) = \frac{\mu (\sigma + \nu)(\mu\nu + 1)}{\nu (1 - \sigma)},}
#' which is finite for \eqn{\sigma < 1} and increases without bound as
#' \eqn{\sigma \to 1}. Both \eqn{\sigma} and \eqn{\nu} add dispersion beyond the
#' \link[=nbinom2]{negative binomial}, which is recovered as \eqn{\sigma \to 0}.
#'
#' Because the mean is pinned to a single parameter, this parameterisation is
#' the more stable of the two to estimate in; see \code{\link{bnbinom}} for the
#' identifiability problem it avoids.
#'
#' There is no distribution function, since the beta-negative binomial
#' distribution function has no closed form and the support is unbounded.
#' One-step-ahead residuals are therefore not available.
#'
#' @references
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [bnbinom], [nbinom2], [betabinom]
#'
#' @param x vector of non-negative counts.
#' @param n number of random values to return (for \code{rbnbinom2}).
#' @param mu mean parameter, must be positive.
#' @param sigma dispersion parameter, must be positive. The variance is finite only for \code{sigma < 1}.
#' @param nu dispersion parameter, must be positive.
#' @param log logical; if \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return
#' \code{dbnbinom2} gives the probability mass function and \code{rbnbinom2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rbnbinom2(5, mu = 4, sigma = 0.4, nu = 0.5)
#' d <- dbnbinom2(x, mu = 4, sigma = 0.4, nu = 0.5)
#' @name bnbinom2
NULL

#' @rdname bnbinom2
#' @export
dbnbinom2 <- function(x, mu, sigma, nu, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(mu <= 0)) stop("mu must be positive.")
    if (any(sigma <= 0)) stop("sigma must be positive.")
    if (any(nu <= 0)) stop("nu must be positive.")
  }

  # potentially escape to RNG; there is no CDF, so OSA is not available
  if (inherits(x, "simref")) {
    return(dGenericSim("dbnbinom2", x = x, mu = mu, sigma = sigma, nu = nu, log = log))
  }
  if (inherits(x, "osa")) {
    stop("Beta-negative binomial does not support OSA residuals.")
  }

  dbnbinom(x, size = 1 / nu, shape1 = 1 / sigma + 1, shape2 = mu * nu / sigma, log = log)
}

#' @rdname bnbinom2
#' @export
rbnbinom2 <- function(n, mu, sigma, nu) {

  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 0)) stop("sigma must be positive.")
  if (any(nu <= 0)) stop("nu must be positive.")

  rbnbinom(n, size = 1 / nu, shape1 = 1 / sigma + 1, shape2 = mu * nu / sigma)
}
