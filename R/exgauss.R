#' Exponentially modified Gaussian distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the exponentially modified Gaussian distribution.
#'
#' @details
#' This implementation of \code{dexgauss} and \code{pexgauss} allows for automatic differentiation with \code{RTMB}.
#' \code{qexgauss} and \code{rexgauss} are reparameterised imports from \code{gamlss.dist::\link[gamlss.dist]{exGAUS}}.
#'
#' If \eqn{X \sim N(\mu, \sigma^2)} and \eqn{Y \sim \text{Exp}(\lambda)}, then
#' \eqn{Z = X + Y} follows the exponentially modified Gaussian distribution with parameters \eqn{\mu}, \eqn{\sigma}, and \eqn{\lambda}.
#'
#' @references
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu mean parameter of the Gaussian part
#' @param sigma standard deviation parameter of the Gaussian part, must be positive.
#' @param lambda rate parameter of the exponential part, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dexgauss} gives the density, \code{pexgauss} gives the distribution function, \code{qexgauss} gives the quantile function, and \code{rexgauss} generates random deviates.
#'
#' @examples
#' x <- rexgauss(1, 1, 2, 2)
#' d <- dexgauss(x, 1, 2, 2)
#' p <- pexgauss(x, 1, 2, 2)
#' q <- qexgauss(p, 1, 2, 2)
#' @name exgauss
NULL

#' @rdname exgauss
#' @export
#' @importFrom RTMB dnorm pnorm
dexgauss <- function(x, mu = 0, sigma = 1, lambda = 1, log = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/exGAUS.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure sigma > 0, lambda > 0
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dexgauss", x=x, mu=mu, sigma=sigma, lambda=lambda, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dexgauss", x=x, mu=mu, sigma=sigma, lambda=lambda, log=log))
  }

  ly <- length(x)
  mu <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)
  lambda <- rep(lambda, length = ly)
  nu <- 1 / lambda

  z <- x - mu - ((sigma * sigma) / nu)

  nu_gr <- greater(nu, 0.05 * sigma) # numerical stability

  logdens <- nu_gr * as.finite.neg(- log(nu) - (z + ((sigma * sigma) / (2 * nu))) / nu + log(1e-300 + RTMB::pnorm(z / sigma))) +
    (1-nu_gr) * RTMB::dnorm(x, mean = mu, sd = sigma, log = TRUE)

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname exgauss
#' @export
#' @importFrom RTMB pnorm
pexgauss <- function(q, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/exGAUS.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    # ensure sigma > 0, lambda > 0
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda must be > 0")
  }

  ly <- length(q)
  mu <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)
  lamba <- rep(lambda, length = ly)
  nu <- 1 / lambda

  z <- q - mu - (sigma^2 / nu)

  nu_gr <- greater(nu, 0.05 * sigma) # numerical stability

  p <- nu_gr * RTMB::pnorm((q-mu)/sigma) - RTMB::pnorm(z/sigma) * exp(((mu+(sigma^2/nu))^2-(mu^2)-2*q*((sigma^2)/nu))/(2*sigma^2)) +
    (1 - nu_gr) * RTMB::pnorm(q, mean = mu, sd = sigma)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname exgauss
#' @export
#' @importFrom gamlss.dist qexGAUS
qexgauss <- function(p, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE) {
  if (!ad_context()) {
    # ensure sigma > 0, lambda > 0
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda must be > 0")
  }

  nu <- 1 / lambda
  gamlss.dist::qexGAUS(p, mu = mu, sigma = sigma, nu = nu,
                       lower.tail = lower.tail, log.p = log.p)
}
#' @rdname exgauss
#' @export
#' @importFrom stats rnorm rexp
rexgauss <- function(n, mu = 0, sigma = 1, lambda = 1) {
  # ensure sigma > 0, lambda > 0
  if (any(sigma <= 0)) stop("sigma must be > 0")
  if (any(lambda <= 0)) stop("lambda must be > 0")

  # Generate n random values from the exponentially modified Gaussian distribution
  rnorm(n, mu, sigma) + rexp(n, rate = lambda)
}
