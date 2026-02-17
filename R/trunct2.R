#' Truncated t distribution with location and scale
#'
#' Density, distribution function, quantile function, and random generation for
#' the truncated t distribution with location \code{mu} and scale \code{sigma}.
#'
#' @details
#' This implementation of \code{dtrunct2} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param df degrees of freedom parameter, must be positive.
#' @param mu location parameter.
#' @param sigma scale parameter, must be positive.
#' @param min,max truncation bounds.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dtrunct2} gives the density, \code{ptrunct2} gives the distribution function,
#' \code{qtrunct2} gives the quantile function, and \code{rtrunct2} generates random deviates.
#'
#' @examples
#' x <- rtrunct2(1, df = 5, mu = 2, sigma = 3, min = -1, max = 5)
#' d <- dtrunct2(x, df = 5, mu = 2, sigma = 3, min = -1, max = 5)
#' p <- ptrunct2(x, df = 5, mu = 2, sigma = 3, min = -1, max = 5)
#' q <- qtrunct2(p, df = 5, mu = 2, sigma = 3, min = -1, max = 5)
#' @name trunct2
NULL

#' @rdname trunct2
#' @export
#' @importFrom RTMB dt
dtrunct2 <- function(x, df, mu = 0, sigma = 1, min = -Inf, max = Inf, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (df <= 0) stop("df must be strictly positive.")
    if (sigma <= 0) stop("sigma must be strictly positive.")
    if (min >= max) stop("min must be less than max.")
  }

  if (inherits(x, "simref")) {
    return(dGenericSim("dtrunct2", x = x, df = df, mu = mu, sigma = sigma,
                       min = min, max = max, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dtrunct2", x = x, df = df, mu = mu, sigma = sigma,
                       min = min, max = max, log = log))
  }

  denom <- pt((max - mu) / sigma, df = df) - pt((min - mu) / sigma, df = df)

  inside <- 0.5 * (1 + sign(x - min) * sign(max - x))

  logdens <- log(inside) + RTMB::dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma) - log(1e-300 + denom)

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname trunct2
#' @export
#' @usage
#' ptrunct2(q, df, mu = 0, sigma = 1, min = -Inf, max = Inf,
#'          lower.tail = TRUE, log.p = FALSE)
ptrunct2 <- function(q, df, mu = 0, sigma = 1, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (df <= 0) stop("df must be strictly positive.")
    if (sigma <= 0) stop("sigma must be strictly positive.")
    if (min >= max) stop("min must be less than max.")
  }

  denom <- pt((max - mu) / sigma, df = df) - pt((min - mu) / sigma, df = df)

  s1 <- sign(q - min)
  val <- (pt((q - mu) / sigma, df = df) - pt((min - mu) / sigma, df = df)) / denom
  s2 <- sign(1 - val)

  p <- 0.5 * (1 + s1 * s2) * val + 0.5 * (1 - s2)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)

  p
}

#' @rdname trunct2
#' @export
#' @usage
#' qtrunct2(p, df, mu = 0, sigma = 1, min = -Inf, max = Inf,
#'          lower.tail = TRUE, log.p = FALSE)
#' @importFrom stats qt
qtrunct2 <- function(p, df, mu = 0, sigma = 1, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE) {

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (df <= 0) stop("df must be strictly positive.")
    if (sigma <= 0) stop("sigma must be strictly positive.")
    if (any(p < 0 | p > 1)) stop("Probabilities must be in [0, 1].")
    if (min >= max) stop("min must be less than max.")
  }

  denom <- pt((max - mu) / sigma, df = df) - pt((min - mu) / sigma, df = df)
  p_untrunc <- p * denom + pt((min - mu) / sigma, df = df)

  mu + sigma * stats::qt(p_untrunc, df = df)
}

#' @rdname trunct2
#' @export
#' @importFrom stats runif qt
rtrunct2 <- function(n, df, mu = 0, sigma = 1, min = -Inf, max = Inf) {

  if (!ad_context()) {
    if (df <= 0) stop("df must be strictly positive.")
    if (sigma <= 0) stop("sigma must be strictly positive.")
    if (min >= max) stop("min must be less than max.")
  }

  u <- runif(n)
  left <- pt((min - mu) / sigma, df = df)
  width <- pt((max - mu) / sigma, df = df) - pt((min - mu) / sigma, df = df)

  mu + sigma * stats::qt(left + u * width, df = df)
}
