#' wrapped Cauchy distribution
#'
#' Density and random generation for the wrapped Cauchy distribution.
#'
#' @details
#' This implementation of \code{dwrpcauchy} allows for automatic differentiation with \code{RTMB}.
#' \code{rwrpcauchy} is simply a wrapper for \code{rwrappedcauchy}imported from \code{circular}.
#'
#' @param x vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param rho concentration parameter of the distribution, must be in the interval from 0 to 1.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of random values to return.
#' @param wrap logical; if \code{TRUE}, generated angles are wrapped to the interval from -pi to pi.
#'
#' @return \code{dwrpcauchy} gives the density and \code{rwrpcauchy} generates random deviates.
#'
#' @examples
#' set.seed(1)
#' x <- rwrpcauchy(10, 0, 0.5)
#' d <- dwrpcauchy(x, 0, 0.5)
#' @name wrpcauchy
NULL

#' @rdname wrpcauchy
#' @export
dwrpcauchy <- function(x, mu = 0, rho, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure rho in [0, 1)
    if (any(rho < 0) || any(rho >= 1)) stop("rho must be in the interval [0, 1).")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dwrpcauchy", x = x, mu = mu, rho = rho, log=log))
  }
  if(inherits(x, "osa")) {
    stop("Wrapped cauchy does not support OSA residuals.")
  }
  rho_sq <- rho * rho

  logdens <- - log(2 * pi) +
    log1p(-rho_sq) -
    log1p(rho_sq - 2 * rho * cos(x - mu))

  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
  }
}

#' @rdname wrpcauchy
#' @export
#' @importFrom circular rwrappedcauchy
rwrpcauchy = function(n, mu = 0, rho, wrap = TRUE) {
  suppressWarnings(
    angles <- as.numeric(rwrappedcauchy(n, mu, rho))
  )

  # if generated angels should be wrapped, i.e. mapped to interval [-pi, pi], do so
  if(wrap){
    angles = (angles + pi) %% (2 * pi) - pi
  }
  angles
}
