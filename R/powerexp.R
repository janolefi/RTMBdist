#' Power Exponential distribution (PE and PE2)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Power Exponential distribution (two versions).
#'
#' @details
#' This implementation of the densities and distribution functions allow for automatic differentiation with \code{RTMB} while the other functions are imported from \code{gamlss.dist} package.
#'
#' For \code{powerexp}, \code{mu} is the mean and \code{sigma} is the standard deviation while this does not hold for \code{powerexp2}.
#'
#' See \code{gamlss.dist::\link[gamlss.dist]{PE}} for more details.
#'
#' @references
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter
#' @param sigma scale parameter, must be positive
#' @param nu shape parameter (real)
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}
#'
#' @return
#' \code{dpowerexp} gives the density, \code{ppowerexp} gives the distribution function, \code{qpowerexp} gives the quantile function, and \code{rpowerexp} generates random deviates.
#'
#' @examples
#' # PE
#' x <- rpowerexp(1, mu = 0, sigma = 1, nu = 2)
#' d <- dpowerexp(x, mu = 0, sigma = 1, nu = 2)
#' p <- ppowerexp(x, mu = 0, sigma = 1, nu = 2)
#' q <- qpowerexp(p, mu = 0, sigma = 1, nu = 2)
#'
#' # PE2
#' x <- rpowerexp2(1, mu = 0, sigma = 1, nu = 2)
#' d <- dpowerexp2(x, mu = 0, sigma = 1, nu = 2)
#' p <- ppowerexp2(x, mu = 0, sigma = 1, nu = 2)
#' q <- qpowerexp2(p, mu = 0, sigma = 1, nu = 2)
#' @name powerexp
NULL

#' @rdname powerexp
#' @export
#' @import RTMB
dpowerexp <- function(x, mu = 0, sigma = 1, nu = 2, log = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/PE.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(nu <= 0)) stop("nu must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dpowerexp", x=x, mu=mu, sigma=sigma, nu=nu, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dpowerexp", x=x, mu=mu, sigma=sigma, nu=nu, log=log))
  }

  log.c <- 0.5 * (-(2/nu) * log(2) + lgamma(1/nu) - lgamma(3/nu))
  c <- exp(log.c)
  z <- (x - mu) / sigma

  logdens <- -log(sigma) + log(nu) - log.c - (0.5*(abs(z/c)^nu)) - (1+(1/nu)) * log(2) - lgamma(1/nu)

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname powerexp
#' @export
#' @usage ppowerexp(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom RTMB pgamma
ppowerexp <- function(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/PE.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(nu <= 0)) stop("nu must be > 0")
  }

  log.c <- 0.5 * (-(2/nu) * log(2) + lgamma(1/nu) - lgamma(3/nu))
  c <- exp(log.c)
  z <- (q - mu) / sigma
  s <- 0.5 * ((abs(z/c))^nu)
  cdf <- 0.5 * (1 + RTMB::pgamma(s, shape=1/nu, scale=1) * sign(z))

  large_nu <- greater(nu, 1e4)

  p <- large_nu * (q-(mu-sqrt(3)*sigma))/(sqrt(12)*sigma) +
    (1-large_nu) * cdf

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}

#' @rdname powerexp
#' @export
#' @usage qpowerexp(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist qPE
qpowerexp <- function(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(nu <= 0)) stop("nu must be > 0")
  }

  gamlss.dist::qPE(p, mu = mu, sigma = sigma, nu = nu,
                   lower.tail = lower.tail, log.p = log.p)
}

#' @rdname powerexp
#' @export
#' @importFrom gamlss.dist rPE
rpowerexp <- function(n, mu = 0, sigma = 1, nu = 2) {

  if (any(sigma <= 0)) stop("sigma must be > 0")
  if (any(nu <= 0)) stop("nu must be > 0")

  gamlss.dist::rPE(n, mu = mu, sigma = sigma, nu = nu)
}

#' @rdname powerexp
#' @export
#' @import RTMB
dpowerexp2 <- function(x, mu = 0, sigma = 1, nu = 2, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args)
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(nu <= 0)) stop("nu must be > 0")
  }

  if (inherits(x, "simref"))
    return(dGenericSim("dpowerexp2", x = x, mu = mu, sigma = sigma, nu = nu, log = log))
  if (inherits(x, "osa"))
    return(dGenericOSA("dpowerexp2", x = x, mu = mu, sigma = sigma, nu = nu, log = log))

  # standardized variable
  z <- (x - mu) / sigma

  # log-density
  logdens <- log(nu) - log(2) - log(sigma) - lgamma(1 / nu) - abs(z)^nu

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname powerexp
#' @export
#' @usage ppowerexp2(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom RTMB pgamma
ppowerexp2 <- function(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/PE2.R
  # and modified to allow for automatic differentiaion

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(nu <= 0)) stop("nu must be > 0")
  }

  z <- (q - mu) / sigma
  s <- abs(z)^nu
  cdf <- 0.5 *(1 + RTMB::pgamma(s, shape=1/nu, scale=1) * sign(z))

  large_nu <- greater(nu, 1e4)

  p <- large_nu * (q-(mu-sigma))/(2*sigma) + (1-large_nu) * cdf

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  return(p)
}

#' @rdname powerexp
#' @export
#' @usage qpowerexp2(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist qPE2
qpowerexp2 <- function(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(nu <= 0)) stop("nu must be > 0")
  }

  gamlss.dist::qPE2(p, mu = mu, sigma = sigma, nu = nu,
                    lower.tail = lower.tail, log.p = log.p)
}

#' @rdname powerexp
#' @export
#' @importFrom gamlss.dist rPE2
rpowerexp2 <- function(n, mu = 0, sigma = 1, nu = 2) {

  if (any(sigma <= 0)) stop("sigma must be > 0")
  if (any(nu <= 0)) stop("nu must be > 0")

  gamlss.dist::rPE2(n, mu = mu, sigma = sigma, nu = nu)
}
