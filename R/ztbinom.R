#' Zero-truncated Binomial distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the zero-truncated Binomial distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' By definition, this distribution only has support on the positive integers (1, ..., size).
#' Any zero-truncated distribution is defined as
#' \deqn{P(X=x | X>0) = P(X=x) / (1 - P(X=0)),}
#' where \eqn{P(X=x)} is the probability mass function of the corresponding untruncated distribution.
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param size number of trials
#' @param prob success probability in each trial
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dztbinom} gives the probability mass function, \code{pztbinom} gives the distribution function, and \code{rztbinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rztbinom(1, size = 10, prob = 0.3)
#' d <- dztbinom(x, size = 10, prob = 0.3)
#' p <- pztbinom(x, size = 10, prob = 0.3)
#' @name ztbinom
NULL

#' @rdname ztbinom
#' @export
#' @importFrom RTMB dbinom
dztbinom <- function(x, size, prob, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(size < 0)) stop("size must be >= 0")
    if (any(prob < 0 | prob > 1)) stop("prob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dztbinom", x=x, size=size, prob=prob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dztbinom", x=x, size=size, prob=prob, log=log))
  }

  log_1m_zprob <- log1p(-dbinom(0, size, prob))  # log(1 - P(X=0))
  logdens <- dbinom(x, size, prob, log = TRUE)

  logdens <- logdens - log_1m_zprob + log(ispos_strict(x))

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname ztbinom
#' @importFrom RTMB pbinom dbinom
#' @export
pztbinom <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    q <- floor(q)  # make sure it's integer-valued
    if (any(size < 0)) stop("size must be >= 0")
    if (any(prob < 0 | prob > 1)) stop("prob must be in [0,1]")
  }

  cdf <- pbinom(q, size, prob)
  p0 <- dbinom(0, size, prob)
  p <- pmax.ad(cdf - p0, 0) / (1 - p0)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}

#' @rdname ztbinom
#' @importFrom stats runif qbinom
#' @importFrom RTMB dbinom
#' @export
rztbinom <- function(n, size, prob) {
  if (any(size < 0)) stop("size must be >= 0")
  if (any(prob < 0 | prob > 1)) stop("prob must be in [0,1]")

  u <- runif(n)
  p0 <- dbinom(0, size, prob)
  x <- qbinom(p0 + (1 - p0) * u, size = size, prob = prob)
  x
}

