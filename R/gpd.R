#' Generalised Pareto distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the generalised Pareto distribution (GPD).
#'
#' @details
#' \code{dgpd} and \code{pgpd} allow for automatic differentiation with \code{RTMB}.
#'
#' With \eqn{z = (x - \mu) / \sigma} the survival function is
#' \deqn{1 - F(x;\,\mu,\sigma,\xi) = \begin{cases} (1 + \xi z)^{-1/\xi} & \xi \neq 0, \\ e^{-z} & \xi = 0, \end{cases}}
#' for \eqn{x \ge \mu}, and the density is
#' \deqn{f(x;\,\mu,\sigma,\xi) = \frac{1}{\sigma} (1 + \xi z)^{-1/\xi - 1}.}
#'
#' This is the limiting distribution of exceedances over a high threshold, so
#' \eqn{\mu} is usually a fixed threshold rather than an estimated parameter.
#' The support is \eqn{x \ge \mu} for \eqn{\xi \ge 0} and
#' \eqn{\mu \le x \le \mu - \sigma/\xi} for \eqn{\xi < 0}. At \eqn{\xi = 0} the
#' distribution is exponential with rate \eqn{1/\sigma}, and at \eqn{\xi > 0}
#' with \eqn{\mu = \sigma/\xi} it is the \link[=pareto]{Pareto} distribution
#' with \eqn{\mu = 1/\xi}.
#'
#' The three cases are covered by one expression, so no branch on the sign or
#' the value of \eqn{\xi} is needed. In particular the derivative with respect
#' to \eqn{\xi} is exact at \eqn{\xi = 0}, which is the usual starting value
#' when the shape is estimated.
#'
#' The threshold itself belongs to the support: \code{dgpd(mu, mu, sigma, xi)}
#' is \eqn{1/\sigma}, as \code{stats::dexp(0, rate)} is \code{rate}. The
#' \code{VGAM}, \code{evd} and \code{extraDistr} implementations return zero
#' there instead.
#'
#' @references
#' Coles, S. (2001) An Introduction to Statistical Modeling of Extreme Values,
#' Springer, doi:10.1007/978-1-4471-3675-0.
#'
#' Pickands, J. (1975) Statistical inference using extreme order statistics.
#' The Annals of Statistics, 3, 119-131.
#'
#' @seealso [gev], [pareto], [frechet]
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter, the threshold below which the density is zero.
#' @param sigma scale parameter, must be positive.
#' @param xi shape parameter (real).
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dgpd} gives the density, \code{pgpd} gives the distribution function, \code{qgpd} gives the quantile function, and \code{rgpd} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rgpd(5, mu = 0, sigma = 1, xi = 0.3)
#' d <- dgpd(x, mu = 0, sigma = 1, xi = 0.3)
#' p <- pgpd(x, mu = 0, sigma = 1, xi = 0.3)
#' q <- qgpd(p, mu = 0, sigma = 1, xi = 0.3)
#' @name gpd
NULL

#' @rdname gpd
#' @export
dgpd <- function(x, mu = 0, sigma = 1, xi = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if(any(sigma <= 0)) stop("sigma must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dgpd", x=x, mu=mu, sigma=sigma, xi=xi, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dgpd", x=x, mu=mu, sigma=sigma, xi=xi, log=log))
  }

  z <- (x - mu) / sigma
  above <- 1 - smaller(z, 0) # 1 for x >= mu, and the threshold itself is included
  insup <- above * greater(1 + xi * z, 0) # xi < 0 also bounds the support above

  zc <- insup * z
  logsurv <- -zc * log1p_over_x(xi * zc) # = log(1 - F), and -z at xi = 0

  logdens <- log(insup) - log(sigma) + (xi + 1) * logsurv

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname gpd
#' @export
pgpd <- function(q, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if(any(sigma <= 0)) stop("sigma must be > 0")
  }

  z <- (q - mu) / sigma
  above <- 1 - smaller(z, 0)
  insup <- above * greater(1 + xi * z, 0)

  zc <- insup * z
  logsurv <- -zc * log1p_over_x(xi * zc)

  # above - insup is 1 only beyond the upper end point, which exists for xi < 0
  p <- insup * (1 - exp(logsurv)) + (above - insup)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  return(p)
}

#' @rdname gpd
#' @export
qgpd <- function(p, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if(any(sigma <= 0)) stop("sigma must be > 0")
  }

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p

  if(!ad_context()) {
    if(any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  n <- max(lengths(list(p, mu, sigma, xi)))
  p <- rep_len(p, n); mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n); xi <- rep_len(xi, n)

  y <- 1 - p # the survival probability, which the two branches invert
  ifelse(xi == 0, mu - sigma * log(y), mu + sigma * (y^(-xi) - 1) / xi)
}

#' @rdname gpd
#' @export
#' @importFrom stats runif
rgpd <- function(n, mu = 0, sigma = 1, xi = 0) {

  if(!ad_context()) {
    if(any(sigma <= 0)) stop("sigma must be > 0")
  }

  n <- ceiling(n)
  p <- runif(n)

  qgpd(p, mu = mu, sigma = sigma, xi = xi)
}
