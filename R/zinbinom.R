#' Zero-inflated negative binomial distribution
#'
#' Probability mass function, distribution function, quantile function, and random generation for
#' the zero-inflated negative binomial distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' \deqn{P(X=k;\,r,p,p_0) = p_0\,\mathbf{1}[k=0] + (1-p_0)\,P_{\mathrm{NB}}(k;\,r,p),}
#' where \eqn{p_0} is \code{zeroprob}.
#'
#' @param x,q vector of (non-negative integer) quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param size size parameter, must be positive.
#' @param prob mean parameter, must be positive.
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzinbinom} gives the density, \code{pzinbinom} gives the distribution function, and \code{rzinbinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzinbinom(1, size = 2, prob = 0.5, zeroprob = 0.5)
#' d <- dzinbinom(x, size = 2, prob = 0.5, zeroprob = 0.5)
#' p <- pzinbinom(x, size = 2, prob = 0.5, zeroprob = 0.5)
#' @name zinbinom
NULL
#' @rdname zinbinom
#' @export
#' @importFrom RTMB dnbinom
dzinbinom <- function(x, size, prob, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure size >= 0, prob in (0,1], zeroprob in [0,1]
    if (any(size <= 0)) stop("size must be > 0")
    if (any(prob <= 0 | prob > 1)) stop("prob must be in (0,1]")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dzinbinom", x = x, size=size, prob=prob, zeroprob=zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzinbinom", x = x, size=size, prob=prob, zeroprob=zeroprob, log=log))
  }

  logdens <- RTMB::dnbinom(x, size = size, prob = prob, log = TRUE)
  logdens <- logspace_add(log(zeroprob) + log(iszero(x)), logdens + log1p(-zeroprob))

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname zinbinom
#' @export
pzinbinom <- function(q, size, prob, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure size >= 0, prob in (0,1], zeroprob in [0,1]
    if (any(size <= 0)) stop("size must be > 0")
    if (any(prob <= 0 | prob > 1)) stop("prob must be in (0,1]")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    q <- floor(q)  # make sure it's integer-valued
  }

  # pnbinom gives 0 for q < 0, so no handling of that case necessary
  p <- zeroprob + (1 - zeroprob) * RTMB::pnbinom(q, size = size, prob = prob)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname zinbinom
#' @importFrom stats runif rnbinom
#' @export
rzinbinom <- function(n, size, prob, zeroprob = 0) {
  # ensure size >= 0, prob in (0,1], zeroprob in [0,1]
  if (any(size <= 0)) stop("size must be > 0")
  if (any(prob <= 0 | prob > 1)) stop("prob must be in (0,1]")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  u <- runif(n)
  res <- ifelse(u < zeroprob, 0, rnbinom(n, size=size, prob=prob))
  return(res)
}



