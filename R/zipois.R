#' Zero-inflated Poisson distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the zero-inflated Poisson distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param lambda vector of (non-negative) means
#' @param zeroprob zero-inflation probability between 0 and 1
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzipois} gives the probability mass function, \code{pzipois} gives the distribution function, and \code{rzipois} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzipois(1, 0.5, 1)
#' d <- dzipois(x, 0.5, 1)
#' p <- pzipois(x, 0.5, 1)
#' @name zipois
NULL
#' @rdname zipois
#' @export
#' @importFrom RTMB logspace_add dpois
dzipois <- function(x, lambda, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure lambda >= 0, zeroprob in [0,1]
    if (any(lambda < 0)) stop("lambda must be >= 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dzipois", x = x, lambda = lambda, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzipois", x = x, lambda = lambda, zeroprob = zeroprob, log=log))
  }

  logdens <- RTMB::dpois(x, lambda = lambda, log = TRUE)
  logdens <- logspace_add(log(zeroprob) + log(iszero(x)), logdens + log1p(-zeroprob))

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname zipois
#' @importFrom RTMB ppois
#' @export
pzipois <- function(q, lambda, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure lambda >= 0, zeroprob in [0,1]
    if (any(lambda < 0)) stop("lambda must be >= 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    q <- floor(q)  # make sure it's integer-valued
  }

  p <- zeroprob + (1 - zeroprob) * ppois(q, lambda)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname zipois
#' @importFrom stats runif rpois
#' @export
rzipois <- function(n, lambda, zeroprob = 0) {
  # ensure lambda >= 0, zeroprob in [0,1]
  if (any(lambda < 0)) stop("lambda must be >= 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  u <- runif(n)
  res <- rep(0, n)
  is_zero <- u < zeroprob
  res[!is_zero] <- rpois(sum(!is_zero), lambda)

  return(res)
}
