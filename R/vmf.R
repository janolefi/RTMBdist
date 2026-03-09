#' von Mises-Fisher distribution
#'
#' Density, distribution function, and random generation for the von Mises-Fisher distribution.
#'
#' @details
#' This implementation of \code{dvmf} allows for automatic differentiation with \code{RTMB}. \code{rvmf} is a reparameterised import from \code{movMF::rmovMF}.
#'
#' @param x unit vector or matrix (with each row being a unit vector) of evaluation points
#' @param mu unit mean vector
#' @param kappa non-negative numeric value for the concentration parameter of the distribution.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of random values to return.
#'
#' @return \code{dvmf} gives the density and \code{rvm} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' # single parameter set
#' mu <- rep(1, 3) / sqrt(3)
#' kappa <- 4
#' x <- rvmf(1, mu, kappa)
#' d <- dvmf(x, mu, kappa)
#'
#' # vectorised over parameters
#' mu <- matrix(mu, nrow = 1)
#' mu <- mu[rep(1,10), ]
#' kappa <- rep(kappa, 10)
#' x <- rvmf(10, mu, kappa)
#' d <- dvmf(x, mu, kappa)
#' @name vmf
NULL

#' @rdname vmf
#' @export
#' @importFrom RTMB besselI
dvmf <- function(x, mu, kappa, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure kappa >= 0
    if (any(kappa < 0)) stop("kappa must be non-negative.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    n <- if (is.matrix(x)) nrow(x) else 1
    x[] <- rvmf(n, mu=mu, kappa=kappa)
    return(0)
  }
  if(inherits(x, "osa")) {
    # return(dGenericOSA("dvm", x = x, mu = mu, kappa = kappa, log=log))
    stop("von Mises-Fisher does not support OSA residuals.")
  }

  # if x or mu are vectors, turn into 1 x p matrices
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(dim(mu))) mu <- matrix(mu, nrow = 1)

  # determine dimension of x
  p <- ncol(x)

  # check if mu has the correct dimension
  if(ncol(mu) != p) stop("x and mu must have the same dimension")

  cprod <- rowSums(mu * x) # t(mu) %*% x for each row fast

  # stable calculation of log(besselI(kappa, p/2-1))
  logI <- log(RTMB::besselI(kappa, p / 2 - 1, expon.scaled = TRUE)) + kappa

  logC <- (p / 2 - 1) * log(kappa) - p / 2 * log(2 * pi) - logI

  logdens <- logC + kappa * cprod

  if(log) return(logdens)

  return(exp(logdens))
}
#' @rdname vmf
#' @export
#' @importFrom movMF rmovMF
rvmf <- function(n, mu, kappa) {
  theta <- mu * kappa

  movMF::rmovMF(n, theta, alpha = 1)
}

