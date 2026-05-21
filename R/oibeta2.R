#' Reparameterised one-inflated beta distribution
#'
#' Density, distribution function, and random generation for
#' the one-inflated beta distribution reparameterised in terms of mean and concentration.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return.
#' @param mu mean parameter, must be in the interval from 0 to 1.
#' @param phi concentration parameter, must be positive.
#' @param oneprob one-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{doibeta2} gives the density, \code{poibeta2} gives the distribution function, and \code{roibeta2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- roibeta2(1, 0.6, 2, 0.5)
#' d <- doibeta2(x, 0.6, 2, 0.5)
#' p <- poibeta2(x, 0.6, 2, 0.5)
#' @name oibeta2
NULL
#' @rdname oibeta2
#' @export
doibeta2 <- function(x, mu, phi, oneprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("doibeta2", x=x, mu=mu, phi=phi, oneprob=oneprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("doibeta2", x=x, mu=mu, phi=phi, oneprob=oneprob, log=log))
  }

  # parameter transformation
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  doibeta(x, shape1 = shape1, shape2 = shape2, oneprob = oneprob, log = log)
}
#' @rdname oibeta2
#' @export
poibeta2 <- function(q, mu, phi, oneprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
  }

  # parameter transformation
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  poibeta(q, shape1 = shape1, shape2 = shape2, oneprob = oneprob,
          lower.tail = lower.tail, log.p = log.p)
}
#' @rdname oibeta2
#' @export
roibeta2 <- function(n, mu, phi, oneprob = 0) {

  # ensure mu in [0,1]
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
  # ensure phi > 0
  if (any(phi <= 0)) stop("phi must be strictly positive.")

  # parameter transformation
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  roibeta(n, shape1 = shape1, shape2 = shape2, oneprob = oneprob)
}
