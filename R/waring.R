#' Waring distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the Waring distribution.
#'
#' @details
#' \code{dwaring} and \code{pwaring} allow for automatic differentiation with \code{RTMB}.
#' The parameterisation follows the \code{WARING} family of the
#' \code{gamlss.dist} package, in which \eqn{\mu} is exactly the mean.
#'
#' \deqn{P(X = k;\, \mu, \sigma) = \frac{B\bigl(k + \tfrac{\mu}{\sigma},\; \tfrac{1}{\sigma} + 2\bigr)}{B\bigl(\tfrac{\mu}{\sigma},\; \tfrac{1}{\sigma} + 1\bigr)}, \quad k = 0, 1, 2, \ldots}
#'
#' The Waring is the beta-geometric: a geometric distribution whose success
#' probability carries a beta prior. It is the two-parameter long-tailed count
#' law behind accident proneness and repeat-buying models, and generalises the
#' \link[=yules]{Yule-Simon} distribution, which is the case \eqn{\sigma = \mu}
#' shifted to start at one. The variance
#' \deqn{\mathrm{Var}(X) = \frac{\mu (\sigma + 1)(\mu + 1)}{1 - \sigma}}
#' is finite only for \eqn{\sigma < 1}, and the tail is a power law throughout.
#'
#' It is exactly the \link[=bnbinom2]{mean-parameterised beta-negative binomial}
#' with \code{nu = 1}. Fixing \code{size} at one is what gives it a closed-form
#' distribution function, which the beta-negative binomial does not have in
#' general, so one-step-ahead residuals are available here but not there.
#'
#' @references
#' Irwin, J. O. (1963) The place of mathematics in medical and biological
#' statistics. Journal of the Royal Statistical Society A, 126, 1-45,
#' doi:10.2307/2982445.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [yules], [bnbinom2], [bnbinom]
#'
#' @param x,q vector of non-negative counts.
#' @param n number of random values to return (for \code{rwaring}).
#' @param mu mean parameter, must be positive.
#' @param sigma dispersion parameter, must be positive. The variance is finite only for \code{sigma < 1}.
#' @param log,log.p logical; if \code{TRUE}, probabilities are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dwaring} gives the probability mass function, \code{pwaring} gives the distribution function, and \code{rwaring} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rwaring(5, mu = 2, sigma = 0.5)
#' d <- dwaring(x, mu = 2, sigma = 0.5)
#' p <- pwaring(x, mu = 2, sigma = 0.5)
#' @name waring
NULL

#' @rdname waring
#' @export
dwaring <- function(x, mu = 2, sigma = 2, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(mu <= 0)) stop("mu must be positive.")
    if (any(sigma <= 0)) stop("sigma must be positive.")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dwaring", x = x, mu = mu, sigma = sigma, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dwaring", x = x, mu = mu, sigma = sigma, log = log))
  }

  dbnbinom2(x, mu = mu, sigma = sigma, nu = 1, log = log)
}

#' @rdname waring
#' @export
pwaring <- function(q, mu = 2, sigma = 2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(mu <= 0)) stop("mu must be positive.")
    if (any(sigma <= 0)) stop("sigma must be positive.")
  }

  shape1 <- 1 / sigma + 1
  shape2 <- mu / sigma

  # the tail sum telescopes: P(X >= m) = B(a, b + m) / B(a, b), which is why
  # this case has a closed form where the beta-negative binomial does not
  qc <- pmax.ad(q, 0) # keep lbeta away from its pole; masked out below

  p <- (1 - exp(lbeta(shape1, shape2 + qc + 1) - lbeta(shape1, shape2))) *
    (1 - smaller(q, 0)) # zero below the support

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}

#' @rdname waring
#' @export
rwaring <- function(n, mu = 2, sigma = 2) {

  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 0)) stop("sigma must be positive.")

  rbnbinom2(n, mu = mu, sigma = sigma, nu = 1)
}
