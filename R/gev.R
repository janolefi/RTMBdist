#' Generalised extreme value distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the generalised extreme value (GEV) distribution.
#'
#' @details
#' \code{dgev} and \code{pgev} allow for automatic differentiation with \code{RTMB}.
#'
#' With \eqn{z = (x - \mu) / \sigma} the distribution function is
#' \deqn{F(x;\,\mu,\sigma,\xi) = \exp\bigl(-t(x)\bigr), \qquad
#'       t(x) = \begin{cases} (1 + \xi z)^{-1/\xi} & \xi \neq 0, \\ e^{-z} & \xi = 0, \end{cases}}
#' and the density is \eqn{f(x) = t(x)^{\xi + 1} e^{-t(x)} / \sigma}.
#'
#' The shape \eqn{\xi} determines the tail and with it the support: for
#' \eqn{\xi > 0} (Frechet case) the distribution is heavy-tailed on
#' \eqn{x > \mu - \sigma/\xi}, for \eqn{\xi < 0} (Weibull case) it is bounded
#' above by \eqn{\mu - \sigma/\xi}, and \eqn{\xi = 0} is the
#' \link[=gumbel]{Gumbel} distribution on the whole real line.
#'
#' The three cases are covered by one expression, so no branch on the sign or
#' the value of \eqn{\xi} is needed. In particular the derivative with respect
#' to \eqn{\xi} is exact at \eqn{\xi = 0}, which is the usual starting value
#' when the shape is estimated.
#'
#' @references
#' Coles, S. (2001) An Introduction to Statistical Modeling of Extreme Values,
#' Springer, doi:10.1007/978-1-4471-3675-0.
#'
#' Jenkinson, A. F. (1955) The frequency distribution of the annual maximum (or
#' minimum) values of meteorological elements. Quarterly Journal of the Royal
#' Meteorological Society, 81, 158-171.
#'
#' @seealso [gumbel], [gpd], [frechet]
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter
#' @param sigma scale parameter, must be positive.
#' @param xi shape parameter (real).
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dgev} gives the density, \code{pgev} gives the distribution function, \code{qgev} gives the quantile function, and \code{rgev} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rgev(5, mu = 0, sigma = 1, xi = 0.2)
#' d <- dgev(x, mu = 0, sigma = 1, xi = 0.2)
#' p <- pgev(x, mu = 0, sigma = 1, xi = 0.2)
#' q <- qgev(p, mu = 0, sigma = 1, xi = 0.2)
#' @name gev
NULL

#' @rdname gev
#' @export
dgev <- function(x, mu = 0, sigma = 1, xi = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if(any(sigma <= 0)) stop("sigma must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dgev", x=x, mu=mu, sigma=sigma, xi=xi, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dgev", x=x, mu=mu, sigma=sigma, xi=xi, log=log))
  }

  z <- (x - mu) / sigma
  insup <- greater(1 + xi * z, 0) # the support is 1 + xi z > 0, and all of R at xi = 0

  # outside the support z is replaced by 0, which leaves log t = 0 rather than
  # an overflowing exponent; the density itself is set to zero by log(insup)
  zc <- insup * z
  logt <- -zc * log1p_over_x(xi * zc) # = -log(1 + xi z) / xi, and -z at xi = 0

  logdens <- log(insup) - log(sigma) + (xi + 1) * logt - exp(logt)

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname gev
#' @export
pgev <- function(q, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if(any(sigma <= 0)) stop("sigma must be > 0")
  }

  z <- (q - mu) / sigma
  insup <- greater(1 + xi * z, 0)

  zc <- insup * z
  logt <- -zc * log1p_over_x(xi * zc)

  # outside the support q lies below the lower end point when xi > 0 and above
  # the upper end point when xi < 0, so the second term is 0 and 1 respectively
  p <- insup * exp(-exp(logt)) + (1 - insup) * smaller(xi, 0)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  return(p)
}

#' @rdname gev
#' @export
qgev <- function(p, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {

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

  y <- -log(p) # = t at the quantile, so the two branches invert t
  ifelse(xi == 0, mu - sigma * log(y), mu + sigma * (y^(-xi) - 1) / xi)
}

#' @rdname gev
#' @export
#' @importFrom stats runif
rgev <- function(n, mu = 0, sigma = 1, xi = 0) {

  if(!ad_context()) {
    if(any(sigma <= 0)) stop("sigma must be > 0")
  }

  n <- ceiling(n)
  p <- runif(n)

  qgev(p, mu = mu, sigma = sigma, xi = xi)
}
