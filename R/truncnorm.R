#' Truncated normal distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the truncated normal distribution.
#'
#' @details
#' This implementation of \code{dtruncnorm} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mean mean parameter, must be positive.
#' @param sd standard deviation parameter, must be positive.
#' @param min,max truncation bounds.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dtruncnorm} gives the density, \code{ptruncnorm} gives the distribution function, \code{qtruncnorm} gives the quantile function, and \code{rtruncnorm} generates random deviates.
#'
#' @examples
#' x <- rtruncnorm(1, mean = 2, sd = 2, min = -1, max = 5)
#' d <- dtruncnorm(x, mean = 2, sd = 2, min = -1, max = 5)
#' p <- ptruncnorm(x, mean = 2, sd = 2, min = -1, max = 5)
#' q <- qtruncnorm(p, mean = 2, sd = 2, min = -1, max = 5)
#' @name truncnorm
NULL

#' @rdname truncnorm
#' @export
#' @importFrom RTMB dnorm pnorm
dtruncnorm <- function(x, mean = 0, sd = 1, min = -Inf, max = Inf, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure sd > 0
    if (sd <= 0) stop("sd must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dtruncnorm", x=x, mean=mean, sd=sd, min=min, max=max, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dtruncnorm", x=x, mean=mean, sd=sd, min=min, max=max, log=log))
  }

  # normalisation constant: probability of being within [min, max]
  denom <- RTMB::pnorm(max, mean, sd) - RTMB::pnorm(min, mean, sd)

  # initialise log-density vector
  logdens <- rep(NaN, length(x))

  # logical vector for values inside the truncation bounds
  # inside <- (x >= min) & (x <= max)
  inside <- 0.5 * (1 + sign(x - min) * sign(max - x))

  logdens <- log(inside) + RTMB::dnorm(x, mean, sd, log = TRUE) - log(1e-300 + denom)

  # compute density only for values inside the bounds
  # logdens[inside] <- RTMB::dnorm(x[inside], mean, sd, log = TRUE) - log(denom)

  # return log-density if requested
  if(log) return(logdens)

  return(exp(logdens))
}

#' @rdname truncnorm
#' @export
#' @usage
#' ptruncnorm(q, mean = 0, sd = 1, min = -Inf, max = Inf,
#'            lower.tail = TRUE, log.p = FALSE)
#' @importFrom RTMB pnorm
ptruncnorm <- function(q, mean = 0, sd = 1, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    # ensure sd > 0
    if (sd <= 0) stop("sd must be strictly positive.")
    # ensure min < max
    if (min >= max) stop("min must be less than max.")
  }

  # normalisation constant: probability of being within [min, max]
  denom <- RTMB::pnorm(max, mean, sd) - RTMB::pnorm(min, mean, sd)

  # Compute standardized CDF
  s1 <- sign(q - min) # for constructing AD-compatible "indicator"
  val <- (RTMB::pnorm(q, mean, sd) - RTMB::pnorm(min, mean, sd)) / denom
  s2 <- sign(1 - val) # for constructing AD-compatible "indicator"

  p <- 0.5 * (1 + s1 * s2) * val + 0.5 * (1 - s2)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)

  return(p)
}

#' @rdname truncnorm
#' @export
#' @usage
#' qtruncnorm(p, mean = 0, sd = 1, min = -Inf, max = Inf,
#'            lower.tail = TRUE, log.p = FALSE)
#' @importFrom RTMB pnorm qnorm
qtruncnorm <- function(p, mean = 0, sd = 1, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE) {
  if (sd <= 0) stop("Standard deviation 'sd' must be positive.")

  # Handle log.p and lower.tail
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    # ensure sd > 0
    if (sd <= 0) stop("sd must be strictly positive.")
    # Check that probabilities are in [0, 1]
    if (any(p < 0 | p > 1)) stop("Probabilities must be between 0 and 1.")
    # ensure min < max
    if (min >= max) stop("min must be less than max.")
  }

  # normalisation constant
  denom <- RTMB::pnorm(max, mean, sd) - RTMB::pnorm(min, mean, sd)

  # Transform p into quantiles of untruncated normal
  p_untrunc <- p * denom + RTMB::pnorm(min, mean, sd)

  # Invert the untruncated normal CDF
  q <- stats::qnorm(p_untrunc, mean, sd)

  return(q)
}

#' @rdname truncnorm
#' @export
#' @importFrom RTMB pnorm qnorm
#' @importFrom stats rnorm runif
rtruncnorm <- function(n, mean = 0, sd = 1, min = -Inf, max = Inf) {

  if (!ad_context()) {
    # ensure sd > 0
    if (sd <= 0) stop("sd must be strictly positive.")
    # ensure min < max
    if (min >= max) stop("min must be less than max.")
  }

  u <- runif(n)
  left <- pnorm((min - mean) / sd)
  right <- pnorm((max - mean) / sd) - pnorm((min - mean) / sd)
  x <- qnorm(left + u * right) * sd + mean
  return(x)
}


