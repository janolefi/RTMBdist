#' Zero-inflated and reparameterised negative binomial distribution
#'
#' Probability mass function, distribution function, quantile function and random generation for
#' the zero-inflated negative binomial distribution reparameterised in terms of mean and size.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' Uses the same density as \code{zinbinom} with \eqn{p = r/(r+\mu)}:
#' \deqn{P(X=k;\,\mu,r,p_0) = p_0\,\mathbf{1}[k=0] + (1-p_0)\,P_{\mathrm{NB}}\!\left(k;\,r,\tfrac{r}{r+\mu}\right).}
#'
#' @param x,q vector of (non-negative integer) quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mu mean parameter, must be positive.
#' @param size size parameter, must be positive.
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzinbinom2} gives the density, \code{pzinbinom2} gives the distribution function, and \code{rzinbinom2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzinbinom2(1, 2, 1, zeroprob = 0.5)
#' d <- dzinbinom2(x, 2, 1, zeroprob = 0.5)
#' p <- pzinbinom2(x, 2, 1, zeroprob = 0.5)
#' @name zinbinom2
NULL
#' @rdname zinbinom2
#' @export
dzinbinom2 <- function(x, mu, size, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure size, mu > 0, zeroprob in [0,1]
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(size <= 0)) stop("size must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dzinbinom2", x = x, mu=mu, size=size, zeroprob=zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzinbinom2", x = x, mu=mu, size=size, zeroprob=zeroprob, log=log))
  }

  # parameter transformation
  prob <- size / (size + mu)

  dzinbinom(x, size = size, prob = prob, zeroprob = zeroprob, log = log)
}
#' @rdname zinbinom2
#' @export
pzinbinom2 <- function(q, mu, size, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {
  if(!ad_context()) {
    # ensure size, mu > 0, zeroprob in [0,1]
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(size <= 0)) stop("size must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    q <- floor(q)  # make sure it's integer-valued
  }

  # parameter transformation
  prob <- size / (size + mu)

  pzinbinom(q, size = size, prob = prob, zeroprob = zeroprob, lower.tail = lower.tail, log.p = log.p)
}
#' @rdname zinbinom2
#' @export
rzinbinom2 <- function(n, mu, size, zeroprob = 0) {
  # ensure mu, size > 0, zeroprob in [0,1]
  if (any(mu <= 0)) stop("mu must be > 0")
  if (any(size <= 0)) stop("size must be > 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  # parameter transformation
  prob <- size / (size + mu)

  rzinbinom(n, size = size, prob = prob, zeroprob = zeroprob)
}

