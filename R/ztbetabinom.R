#' Zero-truncated beta-binomial distribution
#'
#' Probability mass function, distribution function and random generation for the zero-truncated beta-binomial
#' distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' By definition, this distribution only has support on the positive integers (1, ..., n).
#' Any zero-truncated distribution is defined as
#' \deqn{P(X=x | X>0) = P(X=x) / (1 - P(X=0)),}
#' where \eqn{P(X=x)} is the probability mass function of the corresponding
#' untruncated distribution.
#'
#' The distribution function has no closed form and is computed by summing the
#' probability mass function over \eqn{0, \ldots, q}. It is AD-compatible in the
#' parameters, while \code{q} and \code{size} must be numeric data. This is also what
#' one-step-ahead (OSA) residuals via \code{RTMB::\link[RTMB]{oneStepPredict}} need,
#' so these are supported, e.g. with \code{method = "cdf"} and \code{discrete = TRUE}.
#'
#' @seealso [betabinom], [zibetabinom], [hbetabinom], [ztbinom]
#'
#' @param x integer vector of counts
#' @param q vector of quantiles.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le q]}, otherwise \eqn{P[X > q]}.
#' @param log.p logical; if \code{TRUE}, probabilities are returned on the log scale.
#' @param n number of random values to return.
#' @param size number of trials (zero or more)
#' @param shape1,shape2 positive shape parameters of the mixing beta distribution
#' 
#' @param log logical; return log-density if TRUE
#'
#' @return
#' \code{dztbetabinom} gives the probability mass function, \code{pztbetabinom} gives the distribution function, and \code{rztbetabinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rztbetabinom(5, size = 10, shape1 = 2, shape2 = 3)
#' d <- dztbetabinom(x, size = 10, shape1 = 2, shape2 = 3)
#' p <- pztbetabinom(x, size = 10, shape1 = 2, shape2 = 3)
#' @name ztbetabinom
NULL
#' @rdname ztbetabinom
#' @export
dztbetabinom <- function(x, size, shape1, shape2, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(size < 0) || any(size != floor(size))) stop("size must be non-negative integers.")
    if (any(shape1 <= 0) || any(shape2 <= 0)) stop("shape1 and shape2 must be positive.")
  }

  if (inherits(x, "simref")) {
    return(dGenericSim("dztbetabinom", x = x, size = size, shape1 = shape1,
                       shape2 = shape2, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dztbetabinom", x = x, size = size, shape1 = shape1,
                       shape2 = shape2, log = log))
  }

  log_1m_zprob <- log1p(-dbetabinom(0, size, shape1, shape2))  # log(1 - P(X=0))
  # x is clamped because dbetabinom() hits the lgamma pole at negative integers;
  # the ispos_strict() term below removes the clamped values again
  logdens <- dbetabinom(pmax.ad(x, 0), size, shape1, shape2, log = TRUE)

  logdens <- logdens - log_1m_zprob + log(ispos_strict(x))

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname ztbetabinom
#' @export
pztbetabinom <- function(q, size, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(size < 0) || any(size != floor(size))) stop("size must be non-negative integers.")
    if (any(shape1 <= 0) || any(shape2 <= 0)) stop("shape1 and shape2 must be positive.")
  }

  discrete_cdf(dztbetabinom, q, list(size = size, shape1 = shape1, shape2 = shape2),
               upper = size, lower.tail = lower.tail, log.p = log.p)
}
#' @rdname ztbetabinom
#' @export
#' @importFrom stats runif
rztbetabinom <- function(n, size, shape1, shape2) {
  if (any(size < 0) || any(size != floor(size))) stop("size must be non-negative integers.")
  if (any(shape1 <= 0) || any(shape2 <= 0)) stop("shape1 and shape2 must be positive.")

  size <- rep_len(size, n); shape1 <- rep_len(shape1, n); shape2 <- rep_len(shape2, n)

  # rejection sampling: the beta-binomial has no closed-form quantile function
  x <- rbetabinom(n, size, shape1, shape2)
  repeat {
    zero <- which(x == 0)
    if (!length(zero)) break
    x[zero] <- rbetabinom(length(zero), size[zero], shape1[zero], shape2[zero])
  }

  return(x)
}
