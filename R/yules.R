#' Yule-Simon distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the Yule-Simon distribution.
#'
#' @details
#' \code{dyules} and \code{pyules} allow for automatic differentiation with \code{RTMB}.
#'
#' \deqn{P(X = k;\, \rho) = \rho\, B(k,\, \rho + 1), \quad k = 1, 2, \ldots}
#'
#' The Yule-Simon distribution is the classical long-tailed frequency law, used
#' for word counts, city sizes, citation counts and species-per-genus data. Its
#' tail is a power law, \eqn{P(X = k) \sim \rho\,\Gamma(\rho + 1) k^{-(\rho + 1)}},
#' so the mean \eqn{\rho/(\rho - 1)} exists only for \eqn{\rho > 1} and the
#' variance \eqn{\rho^2 / \{(\rho - 1)^2 (\rho - 2)\}} only for \eqn{\rho > 2}.
#'
#' It is the \link[=bnbinom]{beta-negative binomial} with \code{size = 1} and
#' \code{shape2 = 1}, shifted to start at one, and equivalently the
#' \link[=waring]{Waring} distribution with \code{sigma = mu}, shifted the same
#' way. That special case is what gives it a closed-form distribution function,
#' which the beta-negative binomial does not have in general.
#'
#' The support here starts at one, as in \code{VGAM}. The \code{YULE} family of
#' \code{gamlss.dist} instead starts at zero and is parameterised by its mean
#' \eqn{\mu}, which corresponds to \eqn{\rho = (\mu + 1)/\mu}; use
#' \code{dwaring(x, mu, mu)} for that version.
#'
#' @references
#' Simon, H. A. (1955) On a class of skew distribution functions. Biometrika,
#' 42, 425-440, doi:10.1093/biomet/42.3-4.425.
#'
#' @seealso [waring], [bnbinom], [bnbinom2]
#'
#' @param x,q vector of quantiles.
#' @param n number of random values to return (for \code{ryules}).
#' @param shape positive shape parameter \eqn{\rho}.
#' @param log,log.p logical; if \code{TRUE}, probabilities are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dyules} gives the probability mass function, \code{pyules} gives the distribution function, and \code{ryules} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- ryules(5, shape = 2)
#' d <- dyules(x, shape = 2)
#' p <- pyules(x, shape = 2)
#' @name yules
NULL

#' @rdname yules
#' @export
#' @import RTMB
dyules <- function(x, shape = 1, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(shape <= 0)) stop("shape must be positive.")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dyules", x = x, shape = shape, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dyules", x = x, shape = shape, log = log))
  }

  # the beta-negative binomial with size = shape2 = 1, shifted to start at one;
  # that shift also carries the support, since dbnbinom is zero below zero
  dbnbinom(x - 1, size = 1, shape1 = shape, shape2 = 1, log = log)
}

#' @rdname yules
#' @export
pyules <- function(q, shape = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(shape <= 0)) stop("shape must be positive.")
  }

  # the tail sum telescopes: P(X >= m) = B(rho, m) / B(rho, 1), so the
  # distribution function is closed form even though the general
  # beta-negative binomial one is not
  qc <- pmax.ad(q, 0) # keep lbeta away from its pole; masked out below

  p <- (1 - exp(lbeta(shape, qc + 1) - lbeta(shape, 1))) *
    (1 - smaller(q, 1)) # zero below the support

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}

#' @rdname yules
#' @export
ryules <- function(n, shape = 1) {

  if (any(shape <= 0)) stop("shape must be positive.")

  1 + rbnbinom(n, size = 1, shape1 = shape, shape2 = 1)
}
