#' Beta-negative binomial distribution
#'
#' Probability mass function and random generation for the beta-negative
#' binomial distribution.
#'
#' @details
#' \code{dbnbinom} allows for automatic differentiation with \code{RTMB}.
#'
#' The beta-negative binomial arises by giving the success probability of a
#' negative binomial a beta prior, in the same way that the
#' \link[=betabinom]{beta-binomial} does for the binomial:
#' \deqn{P(X = k;\, r, a, b) = \frac{\Gamma(k + r)}{k!\, \Gamma(r)} \frac{B(a + r,\, b + k)}{B(a,\, b)}, \quad k = 0, 1, 2, \ldots}
#'
#' The extra beta layer gives a much heavier tail than the negative binomial:
#' the mean \eqn{rb / (a - 1)} exists only for \eqn{a > 1} and the variance
#' \deqn{\frac{r b (r + a - 1)(b + a - 1)}{(a - 2)(a - 1)^2}}
#' only for \eqn{a > 2}. As \eqn{a \to \infty} with \eqn{b/(a+b)} held fixed the
#' distribution collapses to the negative binomial.
#'
#' Note that the three parameters are only weakly identified: the likelihood has
#' a long ridge along which \code{size} and \code{shape2} trade off against each
#' other, so fitting all three at once needs either a lot of data or a
#' restriction. The \link[=bnbinom2]{mean parameterisation} is usually the more
#' stable one to estimate in.
#'
#' There is no distribution function, since the beta-negative binomial
#' distribution function has no closed form and the support is unbounded.
#' One-step-ahead residuals are therefore not available.
#'
#' @references
#' Johnson, N. L., Kemp, A. W. and Kotz, S. (2005) Univariate Discrete
#' Distributions, 3rd edition, Wiley, doi:10.1002/0471715816.
#'
#' @seealso [bnbinom2], [betabinom], [nbinom2]
#'
#' @param x vector of non-negative counts.
#' @param n number of random values to return (for \code{rbnbinom}).
#' @param size positive number of successes (need not be an integer).
#' @param shape1 positive shape parameter 1 of the beta prior.
#' @param shape2 positive shape parameter 2 of the beta prior.
#' @param log logical; if \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return
#' \code{dbnbinom} gives the probability mass function and \code{rbnbinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rbnbinom(5, size = 3, shape1 = 4, shape2 = 2)
#' d <- dbnbinom(x, size = 3, shape1 = 4, shape2 = 2)
#' @name bnbinom
NULL

#' @rdname bnbinom
#' @export
#' @import RTMB
dbnbinom <- function(x, size, shape1, shape2, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(size <= 0)) stop("size must be positive.")
    if (any(shape1 <= 0) || any(shape2 <= 0)) stop("shape1 and shape2 must be positive.")
  }

  # potentially escape to RNG; there is no CDF, so OSA is not available
  if (inherits(x, "simref")) {
    return(dGenericSim("dbnbinom", x = x, size = size, shape1 = shape1, shape2 = shape2, log = log))
  }
  if (inherits(x, "osa")) {
    stop("Beta-negative binomial does not support OSA residuals.")
  }

  # clamped below zero so that lgamma stays away from its poles; the density
  # itself is set to zero there by the indicator
  xc <- pmax.ad(x, 0)

  logdens <- lgamma(xc + size) - lgamma(xc + 1) - lgamma(size) +
    lbeta(shape1 + size, shape2 + xc) - lbeta(shape1, shape2) +
    log(1 - smaller(x, 0)) # zero below the support

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname bnbinom
#' @export
rbnbinom <- function(n, size, shape1, shape2) {

  if (any(size <= 0)) stop("size must be positive.")
  if (any(shape1 <= 0) || any(shape2 <= 0)) stop("shape1 and shape2 must be positive.")

  n <- ceiling(n)

  # draw the success probability from its beta prior, then the count given it
  prob <- stats::rbeta(n, shape1, shape2)

  stats::rnbinom(n, size = size, prob = prob)
}
