#' Box–Cox t distribution (BCT)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Box–Cox t distribution.
#'
#' @details
#' \code{dbct} and \code{pbct} allow for automatic differentiation with \code{RTMB}.
#' The parameterisation follows the \code{BCT} family of the \code{gamlss.dist} package.
#'
#' The density is
#' \deqn{f(x; \mu, \sigma, \nu, \tau) = \frac{x^{\nu-1}}{\mu^{\nu} \sigma} \frac{f_t(z;\tau)}{F_t\!\left(1/(\sigma|\nu|);\tau\right)}, \quad x > 0,}
#' where \eqn{z = [(x/\mu)^\nu - 1]/(\nu\sigma)} for \eqn{\nu \neq 0} and \eqn{z = \log(x/\mu)/\sigma} for \eqn{\nu = 0}, and \eqn{f_t(\cdot;\tau)} and \eqn{F_t(\cdot;\tau)} are the PDF and CDF of Student's \eqn{t} distribution with \eqn{\tau} degrees of freedom.
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2006) Using the Box-Cox t distribution in GAMLSS to model skewness and kurtosis.
#' Statistical Modelling, 6(3), 209. doi:10.1191/1471082X06st122oa
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [bccg], [bcpe], [skewt]
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
pbct <- function(q, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCT.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(tau <= 0)) stop("tau must be > 0")
  }

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
qbct <- function(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCT.R

  if (!ad_context()) {
    if (any(mu <= 0) || any(sigma <= 0) || any(tau <= 0)) stop("mu, sigma, tau must be > 0")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  # see qbccg: gamlss.dist's ifelse() branching truncates when nu is shorter than p
  n <- max(lengths(list(p, mu, sigma, nu, tau)))
  p <- rep_len(p, n); mu <- rep_len(mu, n); sigma <- rep_len(sigma, n)
  nu <- rep_len(nu, n); tau <- rep_len(tau, n)

  # stats::pt / stats::qt are used explicitly: pt() is masked inside this
  # package by the AD-compatible approximation in R/t2.R
  Fz <- stats::pt(1 / (sigma * abs(nu)), df = tau)
  z <- ifelse(nu <= 0, stats::qt(p * Fz, df = tau), stats::qt(1 - (1 - p) * Fz, df = tau))

  return(inv_boxcox(mu, sigma, nu, z))
}

#' @rdname bct
#' @export
#' @importFrom stats runif
rbct <- function(n, mu = 5, sigma = 0.1, nu = 1, tau = 2) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/BCT.R

  if (!ad_context()) {
    if (length(n) != 1 || !is.finite(n) || n < 0) stop("n must be a non-negative scalar")
    if (any(mu <= 0) || any(sigma <= 0) || any(tau <= 0)) stop("mu, sigma, tau must be > 0")
  }

  n <- ceiling(n)
  p <- runif(n)

  qbct(p, mu = mu, sigma = sigma, nu = nu, tau = tau)
}
