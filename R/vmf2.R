#' Reparameterised von Mises-Fisher distribution
#'
#' Density, distribution function, and random generation for the von Mises-Fisher distribution.
#'
#' @details
#' In this parameterisation, \eqn{\theta = \kappa \mu}, where \eqn{\mu} is a unit vector and \eqn{\kappa} is the concentration parameter.
#'
#' \code{dvmf2} allows for automatic differentiation with \code{RTMB}. \code{rvmf2} is imported from \code{movMF::rmovMF}.
#'
#'
#' @param x unit vector or matrix (with each row being a unit vector) of evaluation points
#' @param theta direction and concentration vector. The direction of \code{theta} determines the mean direction on the sphere.
#' The norm of \code{theta} is the concentration parameter of the distribution.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of random values to return.
#'
#' @return \code{dvmf} gives the density and \code{rvm} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' # single parameter set
#' theta <- c(1,2,3)
#' x <- rvmf2(1, theta)
#' d <- dvmf2(x, theta)
#'
#' # vectorised over parameters
#' theta <- matrix(theta, nrow = 1)
#' theta <- theta[rep(1,10), ]
#' x <- rvmf2(10, theta)
#' d <- dvmf2(x, theta)
#' @name vmf2
NULL

#' @rdname vmf2
#' @export
dvmf2 <- function(x, theta, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    n <- if (is.matrix(x)) nrow(x) else 1
    x[] <- rvmf2(n, theta=theta)
    return(0)
  }

  # if x or theta are vectors, turn into 1 x p matrices
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(dim(theta))) theta <- matrix(theta, nrow = 1)

  # determine dimension of x
  p <- ncol(x)

  # check if mu has the correct dimension
  if(ncol(theta) != p) stop("x and theta must have the same dimension")

  kappa <- sqrt(rowSums(theta^2)) # kappa is the norm of theta

  cprod <- rowSums(theta * x) # t(theta) %*% x for each row fast

  # stable calculation of log(besselI(kappa, p/2-1))
  logI <- log(RTMB::besselI(kappa, p / 2 - 1, expon.scaled = TRUE)) + kappa

  logC <- (p / 2 - 1) * log(kappa) - p / 2 * log(2 * pi) - logI

  logdens <- logC + cprod

  if(log) return(logdens)

  return(exp(logdens))
}
#' @rdname vmf2
#' @export
#' @importFrom movMF rmovMF
rvmf2 <- function(n, theta) {
  movMF::rmovMF(n, theta, alpha = 1)
}

