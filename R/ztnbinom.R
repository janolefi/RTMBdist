#' Zero-truncated Negative Binomial distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the zero-truncated Negative Binomial distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' By definition, this distribution only has support on the positive integers (1, 2, ...).
#' Any zero-truncated distribution is defined as
#' \deqn{P(X=x | X>0) = P(X=x) / (1 - P(X=0)),}
#' where \eqn{P(X=x)} is the probability mass function of the corresponding untruncated distribution.
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dztnbinom} gives the probability mass function, \code{pztnbinom} gives the distribution function, and \code{rztnbinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rztnbinom(1, size = 2, prob = 0.5)
#' d <- dztnbinom(x, size = 2, prob = 0.5)
#' p <- pztnbinom(x, size = 2, prob = 0.5)
#' @name ztnbinom
NULL

#' @rdname ztnbinom
#' @export
#' @importFrom RTMB dnbinom
dztnbinom <- function(x, size, prob, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(size <= 0)) stop("size must be > 0")
    if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
  }

  log_1m_zprob <- log1p(-dnbinom(0, size = size, prob = prob))  # log(1 - P(X=0))
  logdens <- dnbinom(x, size = size, prob = prob, log = TRUE)

  logdens <- logdens - log_1m_zprob + log(ispos_strict(x))

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname ztnbinom
#' @importFrom RTMB pbeta dnbinom
#' @export
pztnbinom <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    q <- floor(q)
    if (any(size <= 0)) stop("size must be > 0")
    if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
  }

  # cdf <- pnbinom(q, size = size, prob = prob)
  # q is clamped because pbeta() is NaN for a negative second shape; below the
  # support the pmax.ad() below then returns 0, as it does for the other families
  cdf <- pbeta(prob, size, pmax.ad(q, 0) + 1)
  p0 <- dnbinom(0, size = size, prob = prob)
  p <- pmax.ad(cdf - p0, 0) / (1 - p0)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}

#' @rdname ztnbinom
#' @importFrom stats runif qnbinom
#' @importFrom RTMB dnbinom
#' @export
rztnbinom <- function(n, size, prob) {
  if (any(size <= 0)) stop("size must be > 0")
  if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")

  u <- runif(n)
  q0 <- dnbinom(0, size, prob)
  x <- qnbinom(q0 + (1-q0)*u, size = size, prob = prob)
  x
}
