#' Pareto distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the pareto distribution.
#'
#' @details
#' This implementation of \code{dpareto} and \code{ppareto} allows for automatic differentiation with \code{RTMB} while the other functions are imported from \code{gamlss.dist} package.
#' See \code{gamlss.dist::\link[gamlss.dist]{PARETO}} for more details.
#'
#' \deqn{f(x;\,\mu) = \frac{\mu}{x^{\mu+1}}, \quad x > 1.}
#'
#' @references
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dpareto} gives the density, \code{ppareto} gives the distribution function, \code{qpareto} gives the quantile function, and \code{rpareto} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rpareto(1, mu = 5)
#' d <- dpareto(x, mu = 5)
#' p <- ppareto(x, mu = 5)
#' q <- qpareto(p, mu = 5)
#' @name pareto
NULL

#' @rdname pareto
#' @export
dpareto <- function(x, mu = 1, log = FALSE) {

  # taken https://github.com/gamlss-dev/gamlss.dist/blob/main/R/PARETO.R
  # and modified to allow for automatic differentiation

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if(any(mu <= 0)) stop("mu must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dpareto", x=x, mu=mu, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dpareto", x=x, mu=mu, log=log))
  }

  logdens <- log(mu) - (mu + 1) * log(x) +
    log(greater(x, 1)) # return - Inf for x <= 1

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname pareto
#' @export
ppareto <- function(q, mu = 1, lower.tail = TRUE, log.p = FALSE) {

  # taken https://github.com/gamlss-dev/gamlss.dist/blob/main/R/PARETO.R
  # and modified to allow for automatic differentiation

  if(!ad_context()) {
    if(any(mu <= 0)) stop("mu must be > 0")
  }

  p <- 1 - q^(-mu)
  p <- p * greater(q, 1)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  return(p)
}

#' @rdname pareto
#' @export
#' @importFrom gamlss.dist qPARETO
qpareto <- function(p, mu = 1, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if(any(mu <= 0)) stop("mu must be > 0")
  }

  gamlss.dist::qPARETO(p, mu = mu, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname pareto
#' @export
#' @importFrom gamlss.dist rPARETO
rpareto <- function(n, mu = 1) {

  if(!ad_context()) {
    if(any(mu <= 0)) stop("mu must be > 0")
  }

  gamlss.dist::rPARETO(n, mu = mu)
}


