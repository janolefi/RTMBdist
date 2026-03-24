#' Zero- and one-inflated beta distribution
#'
#' Density, distribution function, and random generation for
#' the zero-one-inflated beta distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return.
#' @param shape1,shape2 non-negative shape parameters of the beta distribution
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param oneprob zero-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzoibeta} gives the density, \code{pzoibeta} gives the distribution function, and \code{rzoibeta} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzoibeta(1, 2, 2, 0.2, 0.3)
#' d <- dzoibeta(x, 2, 2, 0.2, 0.3)
#' p <- pzoibeta(x, 2, 2, 0.2, 0.3)
#' @name zoibeta
NULL
#' @rdname zoibeta
#' @export
#' @importFrom RTMB dbeta
dzoibeta <- function(x, shape1, shape2, zeroprob = 0, oneprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # shapes positive
    if (any(shape1 <= 0)) stop("shape1 must be positive.")
    if (any(shape2 <= 0)) stop("shape2 must be positive.")
    # ensure zeroprob in [0,1]
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    # ensure oneprob in [0,1]
    if (any(oneprob < 0 | oneprob > 1)) stop("oneprob must be in [0,1]")
    # ensure zeroprob + oneprob <= 1
    if (any(zeroprob + oneprob > 1)) stop("zeroprob + oneprob must be <= 1")
    # ensure x in [0,1)
    # if (any(x < 0 | x >= 1)) stop("x must be in the interval [0, 1).")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dzoibeta", x=x, shape1=shape1, shape2=shape2,
                       zeroprob=zeroprob, oneprob=oneprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzoibeta", x=x, shape1=shape1, shape2=shape2,
                       zeroprob=zeroprob, oneprob=oneprob, log=log))
  }

  logdens <- RTMB::dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE)

  # turn + Inf into finite
  logdens <- as.finite(logdens)

  # add zeromass
  logdens <- RTMB::logspace_add(
    log(iszero(x)) + log(zeroprob),
    log(isnonzero(x)) + log1p(-zeroprob-oneprob) + logdens
  )
  # add onemass
  logdens <- RTMB::logspace_add(
    log(iszero(x-1)) + log(oneprob),
    log(isnonzero(x-1)) + logdens
  )

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname zoibeta
#' @export
#' @usage pzoibeta(q, shape1, shape2, zeroprob = 0, oneprob = 0,
#'          lower.tail = TRUE, log.p = FALSE)
#' @importFrom RTMB pbeta
pzoibeta <- function(q, shape1, shape2, zeroprob = 0, oneprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # shapes positive
    if (any(shape1 <= 0)) stop("shape1 must be positive.")
    if (any(shape2 <= 0)) stop("shape2 must be positive.")
    # ensure zeroprob in [0,1]
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    # ensure oneprob in [0,1]
    if (any(oneprob < 0 | oneprob > 1)) stop("oneprob must be in [0,1]")
    # ensure zeroprob + oneprob <= 1
    if (any(zeroprob + oneprob > 1)) stop("zeroprob + oneprob must be <= 1")
    # ensure x in [0,1)
    # if (any(x < 0 | x >= 1)) stop("x must be in the interval [0, 1).")
  }

  p <- iszero(q) * zeroprob +
    (1-iszero(q)) * (zeroprob + (1 - zeroprob - oneprob) * RTMB::pbeta(q, shape1, shape2)) +
    (1-isneg(q-1)) * oneprob

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname zoibeta
#' @export
#' @importFrom stats rbeta
rzoibeta <- function(n, shape1, shape2, zeroprob = 0, oneprob = 0) {

  # shapes positive
  if (any(shape1 <= 0)) stop("shape1 must be positive.")
  if (any(shape2 <= 0)) stop("shape2 must be positive.")
  # ensure zeroprob in [0,1]
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  # ensure oneprob in [0,1]
  if (any(oneprob < 0 | oneprob > 1)) stop("oneprob must be in [0,1]")
  # ensure zeroprob + oneprob <= 1
  if (any(zeroprob + oneprob > 1)) stop("zeroprob + oneprob must be <= 1")

  u <- runif(n)

  res <- rep(0, n)
  is_zero <- u < zeroprob
  is_one  <- u >= zeroprob & u < zeroprob + oneprob
  is_cont <- !is_zero & !is_one
  res[is_one]  <- 1
  res[is_cont] <- rbeta(sum(is_cont), shape1, shape2)

  return(res)
}
