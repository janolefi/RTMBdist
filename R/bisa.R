#' Birnbaum-Saunders distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the Birnbaum-Saunders (fatigue life) distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' The Birnbaum-Saunders distribution with scale \eqn{\beta} and shape
#' \eqn{\alpha} is the distribution of \eqn{X} for which
#' \deqn{\xi(X) = \frac{1}{\alpha}\left(\sqrt{X/\beta} - \sqrt{\beta/X}\right)
#'       = \frac{2}{\alpha}\sinh\left(\tfrac{1}{2}\log(X/\beta)\right)}
#' is standard normal, so that \eqn{F(x) = \Phi(\xi(x))} and
#' \deqn{f(x) = \phi(\xi(x))\, \frac{x + \beta}{2\alpha\sqrt{\beta}\, x^{3/2}}, \quad x > 0.}
#' The scale \eqn{\beta} is the median, and
#' \eqn{E(X) = \beta(1 + \alpha^2/2)}. The quantile function is
#' \eqn{\beta \exp(2\,\mathrm{asinh}(\alpha z / 2))} with \eqn{z = \Phi^{-1}(p)}.
#'
#' @references
#' Birnbaum, Z. W. and Saunders, S. C. (1969) A new family of life distributions.
#' Journal of Applied Probability, 6, 319-327.
#'
#' @seealso [invgauss], [pgweibull], [gompertz]
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param scale scale parameter \eqn{\beta}, the median, must be positive.
#' @param shape shape parameter \eqn{\alpha}, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dbisa} gives the density, \code{pbisa} gives the distribution function,
#' \code{qbisa} gives the quantile function, and \code{rbisa} generates random deviates.
#'
#' @examples
#' set.seed(1)
#' x <- rbisa(5, scale = 2, shape = 0.5)
#' d <- dbisa(x, 2, 0.5)
#' p <- pbisa(x, 2, 0.5)
#' q <- qbisa(p, 2, 0.5)
#' @name bisa
NULL

#' @rdname bisa
#' @export
dbisa <- function(x, scale = 1, shape, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    check_bisa(scale, shape)
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dbisa", x = x, scale = scale, shape = shape, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dbisa", x = x, scale = scale, shape = shape, log = log))
  }

  # x is kept away from zero, where log(x) would be -Inf, and the density is set to
  # zero below the support by the indicator
  xc <- pmax.ad(x, .Machine$double.xmin)
  xi <- 2 * sinh(0.5 * log(xc / scale)) / shape

  logdens <- dnorm(xi, log = TRUE) + log(xc + scale) - 1.5 * log(xc) -
    0.5 * log(scale) - log(2 * shape) +
    log(1 - smaller(x, 0))

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname bisa
#' @export
pbisa <- function(q, scale = 1, shape, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) check_bisa(scale, shape)

  xi <- 2 * sinh(0.5 * log(q * (1 - smaller(q, 0)) / scale)) / shape # negative q as 0

  # the upper tail as Phi(-xi) rather than 1 - p, which keeps it accurate
  if (!lower.tail) xi <- -xi
  pnorm(xi, log.p = log.p)
}

#' @rdname bisa
#' @export
qbisa <- function(p, scale = 1, shape, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) check_bisa(scale, shape)

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  scale * exp(2 * asinh(shape * qnorm(p) / 2))
}

#' @rdname bisa
#' @export
#' @importFrom stats rnorm
rbisa <- function(n, scale = 1, shape) {

  check_bisa(scale, shape)

  scale * exp(2 * asinh(shape * rnorm(n) / 2))
}

check_bisa <- function(scale, shape) {
  if (any(scale <= 0)) stop("scale must be positive.")
  if (any(shape <= 0)) stop("shape must be positive.")
}
