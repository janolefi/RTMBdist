#' Half-Cauchy distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the half-Cauchy distribution.
#'
#' @details
#' \code{dhalfcauchy} and \code{phalfcauchy} allow for automatic differentiation with \code{RTMB}.
#'
#' The half-Cauchy is the distribution of \eqn{|Y|} for \eqn{Y} Cauchy with scale
#' \eqn{\sigma}:
#' \deqn{f(x;\,\sigma) = \frac{2}{\pi\sigma\bigl(1 + (x/\sigma)^2\bigr)}, \quad x \ge 0,}
#' with distribution function \eqn{F(x) = \frac{2}{\pi}\arctan(x/\sigma)}.
#'
#' It is the standard weakly informative prior for the standard deviation of a
#' hierarchical model, recommended by Gelman (2006) in place of the inverse-gamma:
#' the density is flat and non-zero at the origin, so it does not force the
#' variance component away from zero, while the tail is heavy enough to leave
#' large values unpenalised. That makes it a natural companion to models fitted
#' by the Laplace approximation, where the variance components are exactly the
#' parameters at issue.
#'
#' It has no moments of any order. It is the \link[=halft]{half-t} with
#' \code{df = 1}, and the \link[=foldnorm]{folded normal} with \code{mu = 0} is
#' the corresponding half-normal, which the half-t approaches as \code{df} grows.
#'
#' @references
#' Gelman, A. (2006) Prior distributions for variance parameters in hierarchical
#' models. Bayesian Analysis, 1, 515-534, doi:10.1214/06-BA117A.
#'
#' @seealso [halft], [foldnorm], [trunct]
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param sigma scale parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dhalfcauchy} gives the density, \code{phalfcauchy} gives the distribution function, \code{qhalfcauchy} gives the quantile function, and \code{rhalfcauchy} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rhalfcauchy(5, sigma = 2)
#' d <- dhalfcauchy(x, sigma = 2)
#' p <- phalfcauchy(x, sigma = 2)
#' q <- qhalfcauchy(p, sigma = 2)
#' @name halfcauchy
NULL

#' @rdname halfcauchy
#' @export
#' @import RTMB
dhalfcauchy <- function(x, sigma = 1, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(sigma <= 0)) stop("sigma must be > 0")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dhalfcauchy", x = x, sigma = sigma, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dhalfcauchy", x = x, sigma = sigma, log = log))
  }

  # the squared argument is finite below zero as well, so the indicator alone
  # is enough to carry the support; the origin itself belongs to it
  logdens <- log(2 / pi) - log(sigma) - log1p((x / sigma)^2) +
    log(1 - smaller(x, 0))

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname halfcauchy
#' @export
phalfcauchy <- function(q, sigma = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
  }

  p <- (2 / pi) * atan(q / sigma) * (1 - smaller(q, 0))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}

#' @rdname halfcauchy
#' @export
qhalfcauchy <- function(p, sigma = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  n <- max(lengths(list(p, sigma)))
  p <- rep_len(p, n); sigma <- rep_len(sigma, n)

  # tan(pi/2) is merely large in double precision, so the upper end point is
  # returned explicitly rather than as 1.6e16
  ifelse(p == 1, Inf, sigma * tan(pi * p / 2))
}

#' @rdname halfcauchy
#' @export
#' @importFrom stats runif
rhalfcauchy <- function(n, sigma = 1) {

  if (any(sigma <= 0)) stop("sigma must be > 0")

  n <- ceiling(n)
  p <- runif(n)

  qhalfcauchy(p, sigma = sigma)
}
