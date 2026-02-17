#' Box–Cox Cole and Green distribution (BCCG)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Box–Cox Cole and Green distribution.
#'
#' @details
#' This implementation of \code{dbccg} and \code{pbccg} allows for automatic differentiation with \code{RTMB} while the other functions are imported from \code{gamlss.dist} package.
#' See \code{gamlss.dist::\link[gamlss.dist]{BCCG}} for more details.
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
#' \code{dbccg} gives the density, \code{pbccg} gives the distribution function, \code{qbccg} gives the quantile function, and \code{rbccg} generates random deviates.

#'
#' @examples
#' x <- rbccg(5, mu = 10, sigma = 0.2, nu = 0.5)
#' d <- dbccg(x, mu = 10, sigma = 0.2, nu = 0.5)
#' p <- pbccg(x, mu = 10, sigma = 0.2, nu = 0.5)
#' q <- qbccg(p, mu = 10, sigma = 0.2, nu = 0.5)
#' @name bccg
NULL

#' @rdname bccg
#' @export
#' @import RTMB
dbccg <- function(x, mu = 1, sigma = 0.1, nu = 1, log = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCCG.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(sigma <= 0)) stop("sigma must be > 0")
  }

  # allow simulated references / OSA conventions consistent with your BCPE
  if (inherits(x, "simref")) {
    return(dGenericSim("dbccg", x = x, mu = mu, sigma = sigma, nu = nu, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dbccg", x = x, mu = mu, sigma = sigma, nu = nu, log = log))
  }

  ## length of return value
  # n <- max(length(x), length(mu), length(sigma), length(nu))
  # x <- rep_len(x, n)
  # mu <- rep_len(mu, n)
  # sigma <- rep_len(sigma, n)
  # nu <- rep_len(nu, n)

  ## calculating the pdf
  iz <- iszero(nu)

  # preventing problems with nu = 0
  nu <- nu + .Machine$double.xmin

  z <- (1-iz) * ((((x / mu)^nu) - 1) / (nu * sigma)) +
    iz * (log(x / mu) / sigma)

  logdens <- nu * log(x / mu) - log(sigma) - (z * z) / 2 - log(x) -(log(2*pi)) / 2
  logdens <- logdens - log(RTMB::pnorm(1e-300 + 1 / (sigma * abs(nu))))

  logdens <- log(greater(x, 0)) + logdens

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname bccg
#' @export
#' @usage pbccg(q, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE)
#' @importFrom RTMB pnorm
pbccg <- function(q, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCCG.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(sigma <= 0)) stop("sigma must be > 0")
  }

  ## length of return value
  n <- max(length(q), length(mu), length(sigma), length(nu))
  q <- rep_len(q, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  nu <- rep_len(nu, n)

  z <- rep_len(0, n)
  FYy1 <- rep_len(0, n)
  FYy2 <- rep_len(0, n)
  FYy3 <- rep_len(0, n)

  iz <- iszero(nu)

  z <- (1-iz) * (((q / mu)^nu - 1)/(nu * sigma)) +
    iz * (log(q / mu) / sigma)

  FYy1 <- RTMB::pnorm(z)

  gz <- greater(nu, 0)

  FYy2 <- gz * RTMB::pnorm(-1 / (sigma * abs(nu)))
  FYy3 <- RTMB::pnorm(1 / (sigma * abs(nu)))
  p <- (FYy1 - FYy2) / FYy3

  p <- p * greater(q, 0)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}

#' @rdname bccg
#' @export
#' @usage qbccg(p, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist qBCCG
qbccg <- function(p, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(sigma <= 0)) stop("sigma must be > 0")
  }

  gamlss.dist::qBCCG(p, mu = mu, sigma = sigma, nu = nu,
                     lower.tail = lower.tail, log.p = log.p)
}

#' @rdname bccg
#' @export
#' @importFrom gamlss.dist rBCCG
rbccg <- function(n, mu = 1, sigma = 0.1, nu = 1) {

  if (any(mu <= 0)) stop("mu must be > 0")
  if (any(sigma <= 0)) stop("sigma must be > 0")

  gamlss.dist::rBCCG(n, mu = mu, sigma = sigma, nu = nu)
}
