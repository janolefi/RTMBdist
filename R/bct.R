#' Box–Cox t distribution (BCT)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Box–Cox t distribution.
#'
#' @details
#' This implementation of \code{dbct} and \code{pbct} allows for automatic differentiation with \code{RTMB} while the other functions are imported from \code{gamlss.dist} package.
#' See \code{gamlss.dist::\link[gamlss.dist]{BCT}} for more details.
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
#' @param tau degrees of freedom, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dbct} gives the density, \code{pbct} gives the distribution function, \code{qbct} gives the quantile function, and \code{rbct} generates random deviates.
#'
#' @examples
#' x <- rbct(1, mu = 10, sigma = 0.2, nu = 0.5, tau = 4)
#' d <- dbct(x, mu = 10, sigma = 0.2, nu = 0.5, tau = 4)
#' p <- pbct(x, mu = 10, sigma = 0.2, nu = 0.5, tau = 4)
#' q <- qbct(p, mu = 10, sigma = 0.2, nu = 0.5, tau = 4)
#' @name bct
NULL

#' @rdname bct
#' @export
#' @import RTMB
dbct <- function(x, mu = 5, sigma = 0.1, nu = 1, tau = 2, log = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCT.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(tau <= 0)) stop("tau must be > 0")
  }

  if (inherits(x, "simref")) {
    return(dGenericSim("dbct", x = x, mu = mu, sigma = sigma, nu = nu, tau = tau, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dbct", x = x, mu = mu, sigma = sigma, nu = nu, tau = tau, log = log))
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

  # stabilising log if x = 0
  x <- x + .Machine$double.xmin

  # calculating pdf
  iz <- iszero(nu)

  # preventing problems with nu == 0
  nu <- nu + .Machine$double.xmin

  z <- (1 - iz) * (((x / mu)^nu - 1) / (nu * sigma)) +
    iz * (log(x / mu) / sigma)

  logdens <- (nu-1) * log(x) - nu * log(mu) - log(sigma)
  fTz <- lgamma((tau+1) / 2) - lgamma(tau/2) - 0.5 * log(tau) - lgamma(0.5)
  fTz <- fTz - ((tau+1)/2) * log1p((z*z) / tau)

  logdens <- logdens + fTz - log(1e-300 + pt(1 / (sigma * abs(nu)), df = tau))

  large_tau <- greater(tau, 1e6)
  logdens <- large_tau * dbccg(x, mu, sigma, nu, log = TRUE) +
    (1 - large_tau) * logdens

  logdens <- log(greater(x, 0)) + logdens

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname bct
#' @export
#' @usage pbct(q, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist pBCT
pbct <- function(q, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCT.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(tau <= 0)) stop("tau must be > 0")
  }

  ## length of return value
  n <- max(length(q), length(mu), length(sigma), length(nu), length(tau))
  q <- rep_len(q, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  nu <- rep_len(nu, n)
  tau <- rep_len(tau, n)
  z <- rep_len(0, n)
  FYy <-  FYy1 <- FYy2 <- FYy3 <- rep_len(0, n)

  ##  calculate the cdf
  iz <- iszero(nu)
  z <- (1 - iz) * (((q / mu)^nu - 1) / ((nu + .Machine$double.xmin) * sigma)) +
    iz * (log(q / mu) / sigma)

  FYy1 <- pt(z, tau)
  FYy2 <- greater(nu, 0) * pt(-1 / (sigma * abs(nu)), df = tau)
  FYy3 <- pt(1 / (sigma * abs(nu)), df = tau)

  p <- (FYy1 - FYy2) / FYy3
  p <- p * greater(q, 0)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname bct
#' @export
#' @usage qbct(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist qBCT
qbct <- function(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
    if (mu <= 0 || sigma <= 0 || tau <= 0) stop("mu, sigma, tau must be > 0")
  }

  gamlss.dist::qBCT(p, mu = mu, sigma = sigma, nu = nu, tau = tau,
                    lower.tail = lower.tail, log.p = log.p)
}

#' @rdname bct
#' @export
#' @importFrom gamlss.dist rBCT
rbct <- function(n, mu = 5, sigma = 0.1, nu = 1, tau = 2) {

  if (!ad_context()) {
    if (length(n) != 1 || !is.finite(n) || n < 0) stop("n must be a non-negative scalar")
    if (mu <= 0 || sigma <= 0 || tau <= 0) stop("mu, sigma, tau must be > 0")
  }

  gamlss.dist::rBCT(n, mu = mu, sigma = sigma, nu = nu, tau = tau)
}
