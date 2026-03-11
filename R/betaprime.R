#' Beta prime distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the Beta prime distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' If \eqn{X \sim \text{Beta}(\alpha, \beta)}, then \eqn{\frac{X}{1-X} \sim \text{Betaprime}(\alpha, \beta)}
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param shape1,shape2 non-negative shape parameters of the corresponding Beta distribution
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dbetaprime} gives the density, \code{pbetaprime} gives the distribution function, \code{qbetaprime} gives the quantile function, and \code{rbetaprime} generates random deviates.
#'
#' @examples
#' x <- rbetaprime(1, 2, 1)
#' d <- dbetaprime(x, 2, 1)
#' p <- pbetaprime(x, 2, 1)
#' q <- qbetaprime(p, 2, 1)
#' @name betaprime
NULL

#' @rdname betaprime
#' @export
#' @import RTMB
dbetaprime <- function(x, shape1, shape2, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # shapes positive
    if (any(shape1 <= 0)) stop("shape1 must be positive.")
    if (any(shape2 <= 0)) stop("shape2 must be positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dbetaprime", x=x, shape1=shape1, shape2=shape2, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dbetaprime", x=x, shape1=shape1, shape2=shape2, log=log))
  }

  logB <- lbeta(shape1, shape2)
  logdens <- (shape1 - 1) * log(x) - (shape1 + shape2) * log1p(x) - logB

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname betaprime
#' @export
#' @importFrom RTMB pbeta
pbetaprime <- function(q, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # shapes positive
    if (any(shape1 <= 0)) stop("shape1 must be positive.")
    if (any(shape2 <= 0)) stop("shape2 must be positive.")
  }

  pbeta(q / (1+q), shape1, shape2, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname betaprime
#' @export
#' @importFrom RTMB qbeta
qbetaprime <- function(p, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # shapes positive
    if (any(shape1 <= 0)) stop("shape1 must be positive.")
    if (any(shape2 <= 0)) stop("shape2 must be positive.")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  q <- qbeta(p, shape1, shape2)
  return(q / (1-q))
}

#' @rdname betaprime
#' @export
#' @importFrom stats rbeta
rbetaprime <- function(n, shape1, shape2) {
  if (any(shape1 <= 0)) stop("shape1 must be positive.")
  if (any(shape2 <= 0)) stop("shape2 must be positive.")

  x <- rbeta(n, shape1, shape2)
  return(x / (1 - x))
}
