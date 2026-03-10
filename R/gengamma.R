#' Generalised Gamma distribution (GG)
#'
#' Density, distribution function, quantile function, and random generation for
#' the generalised Gamma distribution.
#'
#' @details
#' This implementation of \code{dgengamma}, \code{pgengamma}, and \code{qgengamma} allows for automatic differentiation with \code{RTMB}.
#'
#' @references
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter, must be positive.
#' @param sigma scale parameter, must be positive.
#' @param nu skewness parameter (real).
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dgengamma} gives the density, \code{pgengamma} gives the distribution function, \code{qgengamma} gives the quantile function, and \code{rgengamma} generates random deviates.
#'
#' @examples
#' x <- rgengamma(5, mu = 4, sigma = 0.5, nu = 0.5)
#' d <- dgengamma(x, mu = 4, sigma = 0.5, nu = 0.5)
#' p <- pgengamma(x, mu = 4, sigma = 0.5, nu = 0.5)
#' q <- qgengamma(p, mu = 4, sigma = 0.5, nu = 0.5)
#' @name gengamma
NULL

#' @rdname gengamma
#' @export
#' @import RTMB
dgengamma <- function(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/GG.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(mu <= 0))  stop("mu must be positive")
    if (any(sigma <= 0))  stop("sigma must be positive")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dgengamma", x=x, mu=mu, sigma=sigma, nu=nu, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dgengamma", x=x, mu=mu, sigma=sigma, nu=nu, log=log))
  }

  # preventing problems with x = 0
  x0 <- x
  x <- x + .Machine$double.xmin

  # preventing problems with nu = 0
  # nu <- nu + .Machine$double.xmin

  log_x_mu <- log(x) - log(mu)
  z <- exp(nu * log_x_mu) # avoid power operator

  small_nu <- smaller(abs(nu), 1e-06)

  ldens_large <- dgamma2(z, 1, sigma * abs(nu), log = TRUE) + log(abs(nu)) + log(z) - log(x)
  ldens_small <- -log(x) - 0.5 * log(2*pi) - log(sigma) - (1 / (2*sigma*sigma)) * log_x_mu*log_x_mu

  logdens <- (1-small_nu) * ldens_large + small_nu * ldens_small

  if(!ad_context()) {
    logdens[x0 <= 0] <- -Inf
  }

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname gengamma
#' @export
#' @import RTMB
pgengamma <- function(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(mu <= 0))  stop("mu must be positive")
    if (any(sigma <= 0))  stop("sigma must be positive")
  }

  # preventing problems with nu = 0
  # nu <- nu + .Machine$double.xmin

  z <- exp(nu * log(q / mu)) # avoid power operator

  small_nu <- smaller(abs(nu), 1e-06)

  cdf_large <- pgamma2(z, 1, sigma*abs(nu), lower.tail = TRUE)

  flip <- isneg(nu) # if nu is negative, flip to upper tail
  cdf_large <- (1-flip) * cdf_large + flip * (1-cdf_large)

  cdf_small <- pnorm(log(q), log(mu), sigma)

  p <- (1-small_nu) * cdf_large + small_nu * cdf_small

  if(!ad_context()) {
    p[q <= 0] <- 0
  }

  # p <- p * greater(p, 0)
  p <- pmin.ad(pmax.ad(p,0),1)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}

#' @rdname gengamma
#' @export
#' @import RTMB
qgengamma <- function(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(mu <= 0))  stop("mu must be positive")
    if (any(sigma <= 0))  stop("sigma must be positive")
  }

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p

  # preventing problems with nu = 0
  # nu <- nu + .Machine$double.xmin

  small_nu <- smaller(abs(nu), 1e-06)

  flip <- isneg(nu) # if nu is negative, flip to upper tail
  p_flip <- (1-flip) * p + flip * (1-p)

  z_large <- qgamma2(p_flip, 1, sigma*abs(nu))
  z_small <- qnorm(p, log(mu), sigma)

  y <- (1-small_nu) * (mu * exp((1/nu) * log(z_large))) + small_nu * exp(z_small)

  if(!ad_context()) {
    y[p == 0] <- 0
    y[p == 1] <- Inf
    y[p <  0] <- NaN
    y[p >  1] <- NaN
  }

  return(y)
}

#' @rdname gengamma
#' @export
#' @import RTMB
#' @importFrom stats runif
rgengamma <- function(n, mu = 1, sigma = 0.5, nu = 1) {

  if (any(mu <= 0))  stop("mu must be positive")
  if (any(sigma <= 0))  stop("sigma must be positive")

  n <- ceiling(n)
  p <- runif(n)

  qgengamma(p, mu=mu, sigma=sigma, nu=nu)
}
