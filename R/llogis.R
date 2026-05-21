#' Log-logistic distribution
#'
#' Density, distribution function, quantile function, and random generation for the log-logistic distribution.
#'
#' @details
#' The log-logistic distribution has density
#' \deqn{f(x;\,\alpha,\beta) = \frac{(\beta/\alpha)\,(x/\alpha)^{\beta-1}}{\bigl(1 + (x/\alpha)^\beta\bigr)^2}, \quad x > 0,}
#' where \eqn{\alpha > 0} is the scale parameter and \eqn{\beta > 0} is the shape parameter.
#' The scale parameter equals the median. Larger \eqn{\beta} gives lighter tails; the distribution
#' has a finite mean only when \eqn{\beta > 1} and a finite variance only when \eqn{\beta > 2}.
#'
#' The log-logistic arises naturally as the distribution of \eqn{X = e^Y} where
#' \eqn{Y \sim \mathrm{Logistic}(\log\alpha,\, 1/\beta)}, which yields the
#' numerically convenient log-density
#' \deqn{\log f(x) = \log\beta - \log x + \log f_{\mathrm{logistic}}\!\bigl(\beta\log(x/\alpha)\bigr).}
#'
#' The CDF is
#' \deqn{F(x;\,\alpha,\beta) = \frac{1}{1 + (x/\alpha)^{-\beta}},}
#' and the quantile function is
#' \deqn{Q(p;\,\alpha,\beta) = \alpha \left(\frac{p}{1-p}\right)^{1/\beta}.}
#'
#' \code{dllogis} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles (\eqn{x > 0}).
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param alpha scale parameter (\eqn{\alpha > 0}); equal to the median.
#' @param beta shape parameter (\eqn{\beta > 0}); controls tail heaviness.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]},
#'   otherwise \eqn{P[X > x]}.
#' @param log.p logical; if \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return
#' \code{dllogis} gives the density, \code{pllogis} gives the distribution function,
#' \code{qllogis} gives the quantile function, and \code{rllogis} generates random deviates.
#'
#' @examples
#' x <- rllogis(5, alpha = 1, beta = 2)
#' dllogis(x, alpha = 1, beta = 2)
#' pllogis(c(0.5, 1, 2), alpha = 1, beta = 2)
#' qllogis(c(0.25, 0.5, 0.75), alpha = 1, beta = 2)
#' @name llogis
NULL

#' @rdname llogis
#' @export
#' @importFrom RTMB dlogis plogis
dllogis <- function(x, alpha = 1, beta = 1, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args)
    if (any(alpha <= 0)) stop("alpha must be strictly positive.")
    if (any(beta  <= 0)) stop("beta must be strictly positive.")
  }

  if (inherits(x, "simref")) {
    return(dGenericSim("dllogis", x = x, alpha = alpha, beta = beta, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dllogis", x = x, alpha = alpha, beta = beta, log = log))
  }

  z <- beta * log(x / alpha)
  logdens <- log(beta) - log(x) + RTMB::dlogis(z, log = TRUE)

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname llogis
#' @export
#' @importFrom RTMB plogis
pllogis <- function(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(alpha <= 0)) stop("alpha must be strictly positive.")
    if (any(beta  <= 0)) stop("beta must be strictly positive.")
  }

  p <- RTMB::plogis(beta * log(q / alpha))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname llogis
#' @export
qllogis <- function(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(alpha <= 0)) stop("alpha must be strictly positive.")
    if (any(beta  <= 0)) stop("beta must be strictly positive.")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  alpha * exp(qlogis(p) / beta)
}

#' @rdname llogis
#' @export
#' @importFrom stats rlogis
rllogis <- function(n, alpha = 1, beta = 1) {
  alpha * exp(rlogis(n, location = 0, scale = 1 / beta))
}
