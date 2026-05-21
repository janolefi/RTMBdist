#' Inverse Gaussian distribution
#'
#' Density, distribution function, and random generation for
#' the inverse Gaussian distribution.
#'
#' @details
#' This implementation of \code{dinvgauss} allows for automatic differentiation with \code{RTMB}.
#' \code{qinvgauss} and \code{rinvgauss} are imported from the \code{statmod} package.
#'
#' \deqn{f(x;\,\mu,\lambda) = \sqrt{\frac{\lambda}{2\pi x^3}} \exp\!\left(-\frac{\lambda(x-\mu)^2}{2\mu^2 x}\right), \quad x > 0.}
#'
#' @param x,q vector of quantiles, must be positive.
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mean location parameter
#' @param shape shape parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' @param ... additional parameter passed to \code{statmod::qinvgauss} for numerical evaluation of the quantile function.
#'
#' @return
#' \code{dinvgauss} gives the density, \code{pinvgauss} gives the distribution function, \code{qinvgauss} gives the quantile function, and \code{rinvgauss} generates random deviates.
#'
#' @examples
#' x <- rinvgauss(1, 1, 0.5)
#' d <- dinvgauss(x, 1, 0.5)
#' p <- pinvgauss(x, 1, 0.5)
#' q <- qinvgauss(p, 1, 0.5)
#' @name invgauss
NULL

#' @rdname invgauss
#' @export
dinvgauss <- function(x, mean = 1, shape = 1, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure mean, shape > 0
    if (any(mean <= 0)) stop("mean must be strictly positive.")
    if (any(shape <= 0)) stop("shape must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dinvgauss", x=x, mean=mean, shape=shape, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dinvgauss", x=x, mean=mean, shape=shape, log=log))
  }

  # avoid issues with zero -> if x == 0, returns +/- Inf
  x <- x + .Machine$double.xmin
  x_cen <- x - mean

  logdens <- 0.5 * log(shape) - 0.5 * log(2 * pi) - 1.5 * log(x) -
    (shape * (x_cen * x_cen)) / (2 * (mean * mean) * x)

  if(log) return(logdens)

  return(exp(logdens))
}
#' @rdname invgauss
#' @export
#' @importFrom RTMB pnorm
pinvgauss <- function(q, mean = 1, shape = 1, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure mean, shape > 0
    if (any(mean <= 0)) stop("mean must be strictly positive.")
    if (any(shape <= 0)) stop("shape must be strictly positive.")
  }

  p <- RTMB::pnorm(sqrt(shape / abs(q)) * (q / mean - 1)) +
    exp(2 * shape / mean) * RTMB::pnorm(-sqrt(shape / abs(q)) * (q / mean + 1))

  # s <- sign(q)
  # p <- 0.5 * (1 + s) * s * p
  p <- ispos(q) * p

  if (!lower.tail) p <- 1 - p
  if (log.p) return(log(p))
  return(p)
}
#' @rdname invgauss
#' @export
#' @importFrom statmod qinvgauss
qinvgauss <- function(p, mean = 1, shape = 1, lower.tail=TRUE, log.p = FALSE, ...) {
  if(!ad_context()) {
    # ensure mean, shape > 0
    if (any(mean <= 0)) stop("mean must be strictly positive.")
    if (any(shape <= 0)) stop("shape must be strictly positive.")
  }

  statmod::qinvgauss(p, mean = mean, shape = shape, lower.tail = lower.tail, log.p = log.p, ...)
}
#' @rdname invgauss
#' @export
#' @importFrom statmod rinvgauss
rinvgauss <- function(n, mean = 1, shape = 1) {
  statmod::rinvgauss(n, mean = mean, shape = shape)
}

