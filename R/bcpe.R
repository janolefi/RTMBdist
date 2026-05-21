# inner functions
f.T <- function(t, tau, log = FALSE){
  log.c <- 0.5 * (-(2 / tau) * log(2) + lgamma(1/tau) - lgamma(3/tau))
  c <- exp(log.c)
  logdens <- log(tau) - log.c - (0.5*(abs(t/c)^tau)) - (1+(1/tau)) * log(2) - lgamma(1/tau)
  if(log) return(logdens)
  return(exp(logdens))
}
F.T <- function(t, tau){
  log.c <- 0.5 * (-(2/tau) * log(2) + lgamma(1/tau) - lgamma(3/tau))
  c <- exp(log.c)
  s <- 0.5 * ((abs(t/c))^tau)
  F.s <- RTMB::pgamma(s, shape = 1/tau, scale = 1)
  cdf <- 0.5*(1 + F.s * sign(t))
  cdf
}

#' Box-Cox Power Exponential distribution (BCPE)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Box-Cox Power Exponential distribution.
#'
#' @details
#' This implementation of \code{dbcpe} and \code{pbcpe} allows for automatic differentiation with \code{RTMB} while the other functions are imported from \code{gamlss.dist} package.
#' See \code{gamlss.dist::\link[gamlss.dist]{BCPE}} for more details.
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
#' @param nu vector of \code{nu} parameter values.
#' @param tau vector of \code{tau} parameter values, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dbcpe} gives the density, \code{pbcpe} gives the distribution function, \code{qbcpe} gives the quantile function, and \code{rbcpe} generates random deviates.
#'
#' @examples
#' x <- rbcpe(1, mu = 5, sigma = 0.1, nu = 1, tau = 1)
#' d <- dbcpe(x, mu = 5, sigma = 0.1, nu = 1, tau = 1)
#' p <- pbcpe(x, mu = 5, sigma = 0.1, nu = 1, tau = 1)
#' q <- qbcpe(p, mu = 5, sigma = 0.1, nu = 1, tau = 1)
#' @name bcpe
NULL

#' @rdname bcpe
#' @export
#' @import RTMB
dbcpe <- function(x, mu = 5, sigma = 0.1, nu = 1, tau = 2, log = FALSE) {

  # taken https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCPE.R
  # and modified to allow for automatic differentiation

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(mu <= 0))  stop("mu must be > 0")
    if (any(sigma <= 0))  stop("sigma must be > 0")
    if (any(tau <= 0))  stop("tau must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dbcpe", x=x, mu=mu, sigma=sigma, nu=nu, tau=tau, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dbcpe", x=x, mu=mu, sigma=sigma, nu=nu, tau=tau, log=log))
  }

  # ## length of return value
  # n <- max(length(x), length(mu), length(sigma), length(nu), length(tau))
  # x <- rep_len(x, n)
  # mu <- rep_len(mu, n)
  # sigma <- rep_len(sigma, n)
  # nu <- rep_len(nu, n)
  # tau <- rep_len(tau, n)
  # z <- rep_len(0, n)
  # FYy <- rep_len(0, n)

  iz <- iszero(nu)

  # preventing problems with nu == 0
  nu <- nu + .Machine$double.xmin

  z <- (1-iz) * (((x / mu)^nu - 1) / (nu * sigma)) +
    iz * (log(x / mu) / sigma)

  logfZ <- f.T(z, tau, log=TRUE) - log(F.T(1 / (sigma * abs(nu)), tau))

  logder <- (nu-1) * log(x) - nu * log(mu) - log(sigma)
  logdens <- logder + logfZ

  logdens <- logdens + log(greater(x, 0))

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname bcpe
#' @export
#' @usage pbcpe(q, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)
#' @import RTMB
pbcpe <- function(q, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE) {

  # taken https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCPE.R
  # and modified to allow for automatic differentiation

  if(!ad_context()) {
    if (any(mu <= 0))  stop("mu must be > 0")
    if (any(sigma <= 0))  stop("sigma must be > 0")
    if (any(tau <= 0))  stop("tau must be > 0")
  }

  ##  calculate the cdf
  iz <- iszero(nu)
  z <- (1-iz) * (((q/mu)^nu-1)/(nu*sigma)) +
    iz * (log(q/mu)/sigma)

  FYy1 <- F.T(z, tau)
  FYy2 <- greater(nu, 0) * F.T(-1/(sigma*abs(nu)), tau)
  FYy3 <- F.T(1 / (sigma*abs(nu)), tau)

  p  <- (FYy1 - FYy2) / FYy3
  p <- p * greater(q, 0)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname bcpe
#' @export
#' @usage qbcpe(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist qBCPE
qbcpe <- function(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if (any(mu < 0))  stop("mu must be > 0")
    if (any(sigma < 0))  stop("sigma must be > 0")
    if (any(tau < 0))  stop("tau must be > 0")
  }

  gamlss.dist::qBCPE(p, mu = mu, sigma = sigma, nu = nu, tau = tau,
                     lower.tail = lower.tail, log.p = log.p)
}
#' @rdname bcpe
#' @export
#' @importFrom gamlss.dist rBCPE
rbcpe <- function(n, mu = 5, sigma = 0.1, nu = 1, tau = 2) {

  if (any(mu <= 0))  stop("mu must be > 0")
  if (any(sigma <= 0))  stop("sigma must be > 0")
  if (any(tau <= 0))  stop("tau must be > 0")

  gamlss.dist::rBCPE(n, mu = mu, sigma = sigma, nu = nu, tau = tau)
}
