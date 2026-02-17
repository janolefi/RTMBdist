#' Truncated t distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the truncated t distribution.
#'
#' @details
#' This implementation of \code{dtrunct} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param df degrees of freedom parameter, must be positive.
#' @param min,max truncation bounds.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dtrunct} gives the density, \code{ptrunct} gives the distribution function,
#' \code{qtrunct} gives the quantile function, and \code{rtrunct} generates random deviates.
#'
#' @examples
#' x <- rtrunct(1, df = 5, min = -1, max = 5)
#' d <- dtrunct(x, df = 5, min = -1, max = 5)
#' p <- ptrunct(x, df = 5, min = -1, max = 5)
#' q <- qtrunct(p, df = 5, min = -1, max = 5)
#' @name trunct
NULL

#' @rdname trunct
#' @export
#' @importFrom RTMB dt
dtrunct <- function(x, df, min = -Inf, max = Inf, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (df <= 0) stop("df must be strictly positive.")
    if (min >= max) stop("min must be less than max.")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dtrunct", x = x, df = df, min = min, max = max, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dtrunct", x = x, df = df, min = min, max = max, log = log))
  }

  # normalisation constant
  denom <- pt(max, df = df) - pt(min, df = df)

  # inside indicator
  inside <- 0.5 * (1 + sign(x - min) * sign(max - x))

  # log-density
  logdens <- log(inside) + RTMB::dt(x, df = df, log = TRUE) - log(1e-300 + denom)

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname trunct
#' @export
#' @usage
#' ptrunct(q, df, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE)
ptrunct <- function(q, df, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (df <= 0) stop("df must be strictly positive.")
    if (min >= max) stop("min must be less than max.")
  }

  denom <- pt(max, df = df) - pt(min, df = df)

  s1 <- sign(q - min)
  val <- (pt(q, df = df) - pt(min, df = df)) / denom
  s2 <- sign(1 - val)

  p <- 0.5 * (1 + s1 * s2) * val + 0.5 * (1 - s2)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)

  p
}

#' @rdname trunct
#' @export
#' @usage
#' qtrunct(p, df, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE)
#' @importFrom stats qt
qtrunct <- function(p, df, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE) {

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (df <= 0) stop("df must be strictly positive.")
    if (any(p < 0 | p > 1)) stop("Probabilities must be in [0, 1].")
    if (min >= max) stop("min must be less than max.")
  }

  denom <- pt(max, df = df) - pt(min, df = df)
  p_untrunc <- p * denom + pt(min, df = df)

  stats::qt(p_untrunc, df = df)
}

#' @rdname trunct
#' @export
#' @importFrom stats runif qt
rtrunct <- function(n, df, min = -Inf, max = Inf) {

  if (!ad_context()) {
    if (df <= 0) stop("df must be strictly positive.")
    if (min >= max) stop("min must be less than max.")
  }

  u <- runif(n)
  left <- pt(min, df = df)
  width <- pt(max, df = df) - pt(min, df = df)

  stats::qt(left + u * width, df = df)
}
