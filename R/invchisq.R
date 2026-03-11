#' Inverse Chi-squared distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the inverse Chi-squared distribution.
#'
#' @details
#' If \eqn{X \sim \text{Chisq}(\nu)}, then \eqn{1/X \sim \text{invChisq}(\nu)}.
#'
#' The inverse Chi-squared distribution with \eqn{\nu} degrees of freedom has
#' density
#' \deqn{f(x) = \frac{(\nu/2)^{\nu/2}}{\Gamma(\nu/2)} x^{-(\nu/2+1)} \exp(-\nu/(2x)), \quad x>0.}
#'
#' This implementation of \code{dinvchisq}, \code{pinvchisq}, and \code{qinvchisq} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles, must be positive.
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param df degrees of freedom (\eqn{\nu > 0})
#' @param scale optional positive scale parameter. Default value of \code{1/df} corresponds to standard inverse gamma
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dinvchisq} gives the density, \code{pinvchisq} gives the distribution function,
#' \code{qinvchisq} gives the quantile function, and \code{rinvchisq} generates random deviates.
#'
#' @examples
#' x <- rinvchisq(1, df = 5)
#' d <- dinvchisq(x, df = 5)
#' p <- pinvchisq(x, df = 5)
#' q <- qinvchisq(p, df = 5)
#' @name invchisq
NULL

#' @rdname invchisq
#' @export
#' @import RTMB
dinvchisq <- function(x, df, scale = 1/df, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure df > 0
    if (any(df <= 0)) stop("df must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dinvchisq", x=x, df=df, scale=scale, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dinvchisq", x=x, df=df, scale=scale, log=log))
  }

  nu2 <- df / 2
  logC <- nu2 * log(nu2 * scale) - lgamma(nu2)
  logdens <- - (nu2 + 1) * log(x) - nu2 * scale / x + logC

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname invchisq
#' @export
#' @import RTMB
pinvchisq <- function(q, df, scale = 1/df, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure df > 0
    if (any(df <= 0)) stop("df must be strictly positive.")
  }

  nu2 <- df / 2
  p <- reggamma(nu2, nu2 * scale / q)

  if(!lower.tail) p <- 1-p
  if(log.p) p <- log(p)
  return(p)
}

#' @rdname invchisq
#' @export
#' @import RTMB
qinvchisq <- function(p, df, scale = 1/df, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure df > 0
    if (any(df <= 0)) stop("df must be strictly positive.")
  }

  nu2 <- df / 2
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p

  q <- RTMB::qgamma(1 - p, shape = nu2, scale = 1)
  nu2 * scale / q
}

#' @rdname invchisq
#' @export
#' @importFrom stats rgamma
rinvchisq <- function(n, df, scale = 1/df) {

  # ensure df > 0
  if(any(df <= 0)) stop("df must be positive.")

  nu2 <- df / 2
  y <- stats::rgamma(n, shape = nu2, rate = 1)
  1 / y
}
