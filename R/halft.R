#' Half-t distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the half-t distribution.
#'
#' @details
#' \code{dhalft} and \code{phalft} allow for automatic differentiation with \code{RTMB},
#' with respect to \code{df} as well as \code{sigma}.
#'
#' The half-t is the distribution of \eqn{|Y|} for \eqn{Y = \sigma T} and
#' \eqn{T} Student t on \eqn{\nu} degrees of freedom:
#' \deqn{f(x;\,\nu,\sigma) = \frac{2}{\sigma} f_T(x/\sigma;\, \nu), \quad x \ge 0,}
#' with distribution function \eqn{F(x) = 2 F_T(x/\sigma;\, \nu) - 1}.
#'
#' Together with the \link[=halfcauchy]{half-Cauchy}, which is the case
#' \code{df = 1}, this is the standard weakly informative prior for the standard
#' deviation of a hierarchical model, recommended by Gelman (2006) in place of
#' the inverse-gamma. The degrees of freedom set how heavy the tail is: small
#' \code{df} leaves large values essentially unpenalised, and as \code{df} grows
#' the distribution approaches the half-normal, which is the
#' \link[=foldnorm]{folded normal} with \code{mu = 0}.
#'
#' The mean
#' \deqn{E(X) = \frac{2\sigma\sqrt{\nu}\,\Gamma\bigl(\tfrac{\nu+1}{2}\bigr)}{\sqrt{\pi}\,(\nu - 1)\,\Gamma\bigl(\tfrac{\nu}{2}\bigr)}}
#' exists only for \eqn{\nu > 1} and the variance
#' \eqn{\sigma^2 \nu / (\nu - 2) - E(X)^2} only for \eqn{\nu > 2}.
#'
#' @references
#' Gelman, A. (2006) Prior distributions for variance parameters in hierarchical
#' models. Bayesian Analysis, 1, 515-534, doi:10.1214/06-BA117A.
#'
#' @seealso [halfcauchy], [foldnorm], [trunct], [t2]
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param df degrees of freedom, must be positive.
#' @param sigma scale parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dhalft} gives the density, \code{phalft} gives the distribution function, \code{qhalft} gives the quantile function, and \code{rhalft} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rhalft(5, df = 3, sigma = 2)
#' d <- dhalft(x, df = 3, sigma = 2)
#' p <- phalft(x, df = 3, sigma = 2)
#' q <- qhalft(p, df = 3, sigma = 2)
#' @name halft
NULL

#' @rdname halft
#' @export
#' @import RTMB
dhalft <- function(x, df, sigma = 1, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(df <= 0)) stop("df must be > 0")
    if (any(sigma <= 0)) stop("sigma must be > 0")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dhalft", x = x, df = df, sigma = sigma, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dhalft", x = x, df = df, sigma = sigma, log = log))
  }

  # the t density is symmetric and finite below zero, so the indicator alone
  # is enough to carry the support; the origin itself belongs to it
  logdens <- log(2) - log(sigma) + RTMB::dt(x / sigma, df = df, log = TRUE) +
    log(1 - smaller(x, 0))

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname halft
#' @export
phalft <- function(q, df, sigma = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(df <= 0)) stop("df must be > 0")
    if (any(sigma <= 0)) stop("sigma must be > 0")
  }

  # bare pt() dispatches to the AD version when the arguments are AD variables
  p <- (2 * pt(q / sigma, df = df) - 1) * (1 - smaller(q, 0))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}

#' @rdname halft
#' @export
qhalft <- function(p, df, sigma = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(df <= 0)) stop("df must be > 0")
    if (any(sigma <= 0)) stop("sigma must be > 0")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  sigma * stats::qt((p + 1) / 2, df = df)
}

#' @rdname halft
#' @export
#' @importFrom stats runif
rhalft <- function(n, df, sigma = 1) {

  if (any(df <= 0)) stop("df must be > 0")
  if (any(sigma <= 0)) stop("sigma must be > 0")

  n <- ceiling(n)
  p <- runif(n)

  qhalft(p, df = df, sigma = sigma)
}
