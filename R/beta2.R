#' Reparameterised beta distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the beta distribution reparameterised in terms of mean and concentration.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' Currently, \code{dbeta} masks \code{RTMB::dbeta} because the latter has a numerically unstable gradient.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mu mean parameter, must be in the interval from 0 to 1.
#' @param phi concentration parameter, must be positive.
#' @param shape1,shape2 non-negative parameters
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param eps for internal use only, don't change.
#'
#' @return
#' \code{dbeta2} gives the density, \code{pbeta2} gives the distribution function, \code{qbeta2} gives the quantile function, and \code{rbeta2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rbeta2(1, 0.5, 1)
#' d <- dbeta2(x, 0.5, 1)
#' p <- pbeta2(x, 0.5, 1)
#' q <- qbeta2(p, 0.5, 1)
#' @name beta2
NULL

#' @rdname beta2
#' @export
#' @import RTMB
dbeta <- function(x, shape1, shape2, log = FALSE, eps = 0) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # shapes positive
    if (any(shape1 <= 0)) stop("shape1 must be positive.")
    if (any(shape2 <= 0)) stop("shape2 must be positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dbeta", x=x, shape1=shape1, shape2=shape2, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dbeta", x=x, shape1=shape1, shape2=shape2, log=log))
  }

  logB <- lbeta.ad(shape1, shape2)
  # logB <- RTMB::lbeta(shape1, shape2)
  logdens <- (shape1 - 1) * log(x + eps) +
    (shape2 - 1) * log1p(-x + eps) - logB

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname beta2
#' @export
dbeta2 <- function(x, mu, phi, log = FALSE, eps = 0) {

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
    return(dGenericSim("dbeta2", x=x, mu=mu, phi=phi, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dbeta2", x=x, mu=mu, phi=phi, log=log))
  }

  # parameter transformations
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  dbeta(x, shape1 = shape1, shape2 = shape2, log = log, eps = eps)
}
#' @rdname beta2
#' @export
#' @importFrom RTMB pbeta
pbeta2 <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  if(!ad_context()) {
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
  }

  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  p <- RTMB::pbeta(q, shape1 = shape1, shape2 = shape2)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname beta2
#' @export
#' @importFrom RTMB qbeta
qbeta2 <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  if(!ad_context()) {
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
  }

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p

  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  RTMB::qbeta(p, shape1 = shape1, shape2 = shape2)
}
#' @rdname beta2
#' @export
#' @importFrom stats rbeta
rbeta2 <- function(n, mu, phi) {
  if(!ad_context()) {
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
  }

  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  stats::rbeta(n, shape1 = shape1, shape2 = shape2)
}
