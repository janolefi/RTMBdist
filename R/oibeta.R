#' One-inflated beta distribution
#'
#' Density, distribution function, and random generation for
#' the one-inflated beta distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return.
#' @param shape1,shape2 non-negative shape parameters of the beta distribution
#' @param oneprob one-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{doibeta} gives the density, \code{poibeta} gives the distribution function, and \code{roibeta} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- roibeta(1, 2, 2, 0.5)
#' d <- doibeta(x, 2, 2, 0.5)
#' p <- poibeta(x, 2, 2, 0.5)
#' @name oibeta
NULL
#' @rdname oibeta
#' @export
#' @importFrom RTMB dbeta
doibeta <- function(x, shape1, shape2, oneprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # shapes non-negative
    if (any(shape1 <= 0)) stop("shape1 must be positive.")
    if (any(shape2 <= 0)) stop("shape2 must be positive.")
    # ensure oneprob in [0,1]
    if (any(oneprob < 0 | oneprob > 1)) stop("oneprob must be in [0,1]")
    # ensure x in (0,1]
    # if (any(x <= 0 | x > 1)) stop("x must be in the interval (0, 1].")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("doibeta", x=x, shape1=shape1, shape2=shape2, oneprob=oneprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("doibeta", x=x, shape1=shape1, shape2=shape2, oneprob=oneprob, log=log))
  }

  logdens <- dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE,
                   eps = .Machine$double.xmin)
  logdens <- log_zi(x-1, logdens, oneprob) # use zi function for one inflation by shifting x

  # making sure x == 0 evaluates to -Inf
  logdens <- logdens + log(1-iszero(x))

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname oibeta
#' @export
#' @importFrom RTMB pbeta
poibeta <- function(q, shape1, shape2, oneprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # shapes non-negative
    if (any(shape1 <= 0)) stop("shape1 must be positive.")
    if (any(shape2 <= 0)) stop("shape2 must be positive.")
    # ensure oneprob in [0,1]
    if (any(oneprob < 0 | oneprob > 1)) stop("oneprob must be in [0,1]")
    # ensure x in (0,1]
    # if (any(q <= 0 | q > 1)) stop("q must be in the interval (0, 1].")
  }

  p <- ((1 - oneprob) * RTMB::pbeta(q, shape1, shape2)) +
    ispos_strict(q-1) * oneprob

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname oibeta
#' @export
#' @importFrom stats rbeta
roibeta <- function(n, shape1, shape2, oneprob = 0) {

  # shapes non-negative
  if (any(shape1 <= 0)) stop("shape1 must be positive.")
  if (any(shape2 <= 0)) stop("shape2 must be positive.")
  # ensure oneprob in [0,1]
  if (any(oneprob < 0 | oneprob > 1)) stop("oneprob must be in [0,1]")

  u <- runif(n)
  res <- rep(1, n)
  is_one <- u < oneprob
  res[!is_one] <- rbeta(sum(!is_one), shape1, shape2)

  return(res)
}
