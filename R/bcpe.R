# inner functions
# f.T, F.T and q.T are the density, cdf and quantile function of the standardised
# power exponential, taken from
# https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCPE.R
# (f.T and F.T modified to allow for automatic differentiation)
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
# quantile function of the standardised power exponential; inverse of F.T
q.T <- function(p, tau){
  log.c <- 0.5 * (-(2/tau) * log(2) + lgamma(1/tau) - lgamma(3/tau))
  c <- exp(log.c)
  s <- stats::qgamma((2 * p - 1) * sign(p - 0.5), shape = 1/tau, scale = 1)
  z <- sign(p - 0.5) * ((2 * s)^(1/tau)) * c
  z
}

#' Box-Cox Power Exponential distribution (BCPE)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Box-Cox Power Exponential distribution.
#'
#' @details
#' \code{dbcpe} and \code{pbcpe} allow for automatic differentiation with \code{RTMB}.
#' The parameterisation follows the \code{BCPE} family of the \code{gamlss.dist} package.
#'
#' The density is
#' \deqn{f(x; \mu, \sigma, \nu, \tau) = \frac{x^{\nu-1}}{\mu^{\nu} \sigma} \frac{f_T(z;\tau)}{F_T\!\left(1/(\sigma|\nu|);\tau\right)}, \quad x > 0,}
#' where \eqn{z = [(x/\mu)^\nu - 1]/(\nu\sigma)} for \eqn{\nu \neq 0} and \eqn{z = \log(x/\mu)/\sigma} for \eqn{\nu = 0}, and \eqn{f_T(\cdot;\tau)} and \eqn{F_T(\cdot;\tau)} are the PDF and CDF of the power exponential (PE) distribution with shape \eqn{\tau}.
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2004) Smooth centile curves for skew and kurtotic data modelled using the Box-Cox Power Exponential distribution.
#' Statistics in Medicine, 23, 3053-3076.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [bccg], [bct], [powerexp]
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
qbcpe <- function(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCPE.R

  if(!ad_context()) {
    if (any(mu < 0))  stop("mu must be > 0")
    if (any(sigma < 0))  stop("sigma must be > 0")
    if (any(tau < 0))  stop("tau must be > 0")
  }

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p

  if(!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  # see qbccg: gamlss.dist's ifelse() branching truncates when nu is shorter than p
  n <- max(lengths(list(p, mu, sigma, nu, tau)))
  p <- rep_len(p, n); mu <- rep_len(mu, n); sigma <- rep_len(sigma, n)
  nu <- rep_len(nu, n); tau <- rep_len(tau, n)

  FT <- F.T(1 / (sigma * abs(nu)), tau)
  za <- ifelse(nu < 0, q.T(p * FT, tau), q.T(1 - (1 - p) * FT, tau))
  za <- ifelse(nu == 0, q.T(p, tau), za)

  return(inv_boxcox(mu, sigma, nu, za))
}
#' @rdname bcpe
#' @export
#' @importFrom stats runif
rbcpe <- function(n, mu = 5, sigma = 0.1, nu = 1, tau = 2) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCPE.R

  if (any(mu <= 0))  stop("mu must be > 0")
  if (any(sigma <= 0))  stop("sigma must be > 0")
  if (any(tau <= 0))  stop("tau must be > 0")

  n <- ceiling(n)
  p <- runif(n)

  qbcpe(p, mu = mu, sigma = sigma, nu = nu, tau = tau)
}
