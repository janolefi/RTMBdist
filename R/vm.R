#' von Mises distribution
#'
#' Density, distribution function, and random generation for the von Mises distribution.
#'
#' @details
#' This implementation of \code{dvm} allows for automatic differentiation with \code{RTMB}.
#' \code{rvm} and \code{pvm} are simply wrappers of the corresponding functions from \code{circular}.
#'
#' \deqn{f(x;\,\mu,\kappa) = \frac{\exp(\kappa\cos(x-\mu))}{2\pi\, I_0(\kappa)},}
#' where \eqn{I_0} is the modified Bessel function of the first kind of order 0.
#'
#' @param x,q vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param kappa non-negative numeric value for the concentration parameter of the distribution.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of random values to return.
#' @param tol the precision in evaluating the distribution function
#' @param from value from which the integration for CDF starts. If \code{NULL}, is set to \code{mu - pi}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param log.p logical; if \code{TRUE}, probabilities are returned as \eqn{\log(p)}.
#' @param wrap logical; if \code{TRUE}, generated angles are wrapped to the interval from -pi to pi.
#'
#' @return \code{dvm} gives the density, \code{pvm} gives the distribution function, and \code{rvm} generates random deviates.
#'
#' @examples
#' set.seed(1)
#' x <- rvm(10, 0, 1)
#' d <- dvm(x, 0, 1)
#' p <- pvm(x, 0, 1)
#' @name vm
NULL

#' @rdname vm
#' @export
#' @importFrom RTMB besselI
dvm = function(x, mu = 0, kappa = 1, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure kappa >= 0
    if (any(kappa < 0)) stop("kappa must be non-negative.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dvm", x = x, mu = mu, kappa = kappa, log=log))
  }
  if(inherits(x, "osa")) {
    # return(dGenericOSA("dvm", x = x, mu = mu, kappa = kappa, log=log))
    stop("von Mises does not support OSA residuals.")
  }

  # stable calculation of log(besselI(kappa, 0))
  logI0 <- log(RTMB::besselI(kappa, 0, expon.scaled = TRUE)) + kappa

  logdens <- -log(2 * pi) - logI0 + kappa * cos(x - mu)

  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
  }
}

#' @rdname vm
#' @export
#' @importFrom circular pvonmises
pvm = function(q, mu = 0, kappa = 1, from = NULL, tol = 1e-20,
               lower.tail = TRUE, log.p = FALSE) {
  # NA handling
  ind = which(!is.na(q))

  if(is.matrix(mu)){
    mu = mu[ind,]
  }
  if(is.matrix(kappa)){
    kappa = kappa[ind,]
  }

  probs = numeric(length(q))

  suppressWarnings(
    probs[ind] <- pvonmises(q[ind], mu, kappa, from = from, tol = tol)
  )

  probs[-ind] = NA

  probs <- as.numeric(probs)
  if (!lower.tail) probs <- 1 - probs
  if (log.p) probs <- log(probs)
  probs
}

#' @rdname vm
#' @export
#' @importFrom circular rvonmises
rvm = function(n, mu = 0, kappa = 1, wrap = TRUE) {
  suppressWarnings(
    angles <- as.numeric(rvonmises(n, mu, kappa))
  )

  # if generated angels should be wrapped, i.e. mapped to interval [-pi, pi], do so
  if(wrap){
    angles = (angles + pi) %% (2 * pi) - pi
  }
  angles
}
