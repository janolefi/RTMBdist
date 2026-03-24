#' Zero-inflated binomial distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the zero-inflated binomial distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return.
#' @param size number of trials (zero or more).
#' @param prob probability of success on each trial.
#' @param zeroprob zero-inflation probability between 0 and 1
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzibinom} gives the probability mass function, \code{pzibinom} gives the distribution function, and \code{rzibinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzibinom(1, size = 10, prob = 0.5, zeroprob = 0.5)
#' d <- dzibinom(x, size = 10, prob = 0.5, zeroprob = 0.5)
#' p <- pzibinom(x, size = 10, prob = 0.5, zeroprob = 0.5)
#' @name zibinom
NULL
#' @rdname zibinom
#' @export
#' @importFrom RTMB logspace_add dbinom
dzibinom <- function(x, size, prob, zeroprob = 0, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure size positive integer >= 1, prob in [0,1], zeroprob in [0,1]
    if (any(size < 0 | floor(size) != size)) stop("size must be a non-negative integer")
    if (any(prob < 0 | prob > 1)) stop("prob must be in [0,1]")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dzibinom", x = x, size = size, prob=prob, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzibinom", x = x, size = size, prob=prob, zeroprob = zeroprob, log=log))
  }

  logdens <- RTMB::dbinom(x, size = size, prob = prob, log = TRUE)
  logdens <- log_zi_discrete(x, logdens, zeroprob)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname zibinom
#' @importFrom RTMB pbinom
#' @export
pzibinom <- function(q, size, prob, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    # ensure size positive integer >= 1, prob in [0,1], zeroprob in [0,1]
    if (any(size < 0 | floor(size) != size)) stop("size must be a non-negative integer")
    if (any(prob < 0 | prob > 1)) stop("prob must be in [0,1]")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    q <- floor(q)  # make sure it's integer-valued
  }

  p <- iszero(q) * zeroprob +
    ispos(q) * (zeroprob + (1 - zeroprob) * RTMB::pbinom(q, size=size, prob=prob))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname zibinom
#' @importFrom stats runif rbinom
#' @export
rzibinom <- function(n, size, prob, zeroprob = 0) {
  # ensure size positive integer >= 1, prob in [0,1], zeroprob in [0,1]
  if (any(size < 0 | floor(size) != size)) stop("size must be a non-negative integer")
  if (any(prob < 0 | prob > 1)) stop("prob must be in [0,1]")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  u <- runif(n)
  res <- rep(0, n)
  is_zero <- u < zeroprob
  res[!is_zero] <- rbinom(sum(!is_zero), size=size, prob=prob)

  return(res)
}
