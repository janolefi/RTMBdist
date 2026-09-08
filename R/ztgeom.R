#' Zero-truncated geometric distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the zero-truncated geometric distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' By definition, this distribution only has support on the positive integers (1, 2, ...).
#' Any zero-truncated distribution is defined as
#' \deqn{P(X=x | X>0) = P(X=x) / (1 - P(X=0)),}
#' where \eqn{P(X=x)} is the probability mass function of the corresponding untruncated distribution.
#' For the geometric with success probability \eqn{\pi} this gives
#' \deqn{P(X=x | X>0) = \pi\,(1-\pi)^{x-1}, \quad x = 1, 2, \ldots}
#'
#' @seealso [zigeom], [hgeom], [ztnbinom], [ztpois]
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param prob probability of success in each trial, in (0,1)
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dztgeom} gives the probability mass function, \code{pztgeom} gives the distribution function, and \code{rztgeom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rztgeom(5, prob = 0.3)
#' d <- dztgeom(x, prob = 0.3)
#' p <- pztgeom(x, prob = 0.3)
#' @name ztgeom
NULL
#' @rdname ztgeom
#' @export
#' @import RTMB
dztgeom <- function(x, prob, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dztgeom", x = x, prob = prob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dztgeom", x = x, prob = prob, log=log))
  }

  log_1m_zprob <- log1p(-prob)  # log(1 - P(X = 0)), and P(X = 0) = prob
  logdens <- dgeom(x, prob, log = TRUE)

  logdens <- logdens - log_1m_zprob + log(ispos_strict(x))

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname ztgeom
#' @export
pztgeom <- function(q, prob, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
    q <- floor(q)  # make sure it's integer-valued
  }

  cdf <- pgeom(q, prob)
  p0 <- prob        # P(X = 0)
  p <- pmax.ad(cdf - p0, 0) / (1 - p0)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname ztgeom
#' @export
#' @importFrom stats runif
rztgeom <- function(n, prob) {
  if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")

  u <- runif(n)
  stats::qgeom(prob + (1 - prob) * u, prob)
}
