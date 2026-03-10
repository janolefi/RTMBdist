#' Skellam distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the Skellam distribution.
#'
#' The Skellam distribution is the distribution of the difference of two Poisson random variables. Specifically, if \eqn{X_1 \sim \text{Pois}(\mu_1)} and \eqn{X_2 \sim \text{Pois}(\mu_2)}, then \eqn{X_1 - X_2 \sim \text{Skellam}(\mu_1, \mu_2)}.
#'
#' @details
#' This implementation of \code{dskellam} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x integer vector of counts
#' @param n number of random values to return.
#' @param mu1,mu2 Poisson means
#' @param log logical; return log-density if TRUE
#'
#' @return
#' \code{dskellam} gives the probability mass function and \code{rskellam} generates random deviates.
#'
#' @examples
#' x <- rskellam(1, 2, 3)
#' d <- dskellam(x, 2, 3)
#' @name skellam
NULL
#' @rdname skellam
#' @export
#' @import RTMB
dskellam <- function(x, mu1, mu2, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure lambda, phi > 0
    if (any(mu1 <= 0)) stop("mu1 must be > 0")
    if (any(mu2 <= 0)) stop("mu2 must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dskellam", x = x, mu1 = mu1, mu2 = mu2, log=log))
  }
  if(inherits(x, "osa")) {
    stop("Currently, skellam does not support OSA residuals.")
    # return(dGenericOSA("dskellam", x = x, lambda = lambda, phi = phi, log=log))
  }

  val <- 2 * sqrt(mu1 * mu2)
  logI <- log(besselI(val, x, expon.scaled = TRUE)) + val

  logdens <- -mu1 - mu2 + (x / 2) * (log(mu1) - log(mu2)) + logI

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname skellam
#' @export
rskellam <- function(n, mu1, mu2) {

  if (any(mu1 <= 0)) stop("mu1 must be > 0")
  if (any(mu2 <= 0)) stop("mu2 must be > 0")

  x1 <- rpois(n, mu1)
  x2 <- rpois(n, mu2)

  return(x1 - x2)
}

