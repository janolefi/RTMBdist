#' Zero-inflated log normal distribution
#'
#' Density, distribution function, and random generation for
#' the zero-inflated log normal distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return
#' @param meanlog,sdlog mean and standard deviation of the distribution on the log scale with default values of 0 and 1 respectively.
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzilnorm} gives the density, \code{pzilnorm} gives the distribution function, and \code{rzilnorm} generates random deviates.
#'
#' @examples
#' x <- rzilnorm(1, 1, 1, 0.5)
#' d <- dzilnorm(x, 1, 1, 0.5)
#' p <- pzilnorm(x, 1, 1, 0.5)
#' @name zilnorm
NULL

#' @rdname zilnorm
#' @export
#' @importFrom RTMB dlnorm
dzilnorm = function(x, meanlog = 0, sdlog = 1, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure sdlog > 0, zeroprob in [0,1]
    if (any(sdlog <= 0)) stop("sdlog must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dzilnorm", x=x, meanlog=meanlog, sdlog=sdlog, zeroprob=zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzilnorm", x=x, meanlog=meanlog, sdlog=sdlog, zeroprob=zeroprob, log=log))
  }

  eps <- .Machine$double.xmin # so that gradient is not NaN bc -Inf * 0

  logdens <- RTMB::dlnorm(x + eps, meanlog = meanlog, sdlog = sdlog, log = TRUE)
  logdens <- log_zi(x, logdens, zeroprob)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname zilnorm
#' @export
#' @usage pzilnorm(q, meanlog = 0, sdlog = 1, zeroprob = 0,
#'          lower.tail = TRUE, log.p = FALSE)
pzilnorm <- function(q, meanlog = 0, sdlog = 1, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure sdlog > 0, zeroprob in [0,1]
    if (any(sdlog <= 0)) stop("sdlog must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  eps <- .Machine$double.xmin # so that gradient is not NaN bc -Inf * 0

  p <- iszero(q) * zeroprob +
    ispos_strict(q) * (zeroprob + (1 - zeroprob) * plnorm(q + eps, meanlog = meanlog, sdlog = sdlog))

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname zilnorm
#' @importFrom stats runif rlnorm
#' @export
rzilnorm <- function(n, meanlog = 0, sdlog = 1, zeroprob = 0) {
  # ensure sdlog > 0, zeroprob in [0,1]
  if (any(sdlog <= 0)) stop("sdlog must be > 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  u <- runif(n)
  res <- rep(0, n)
  is_zero <- u < zeroprob
  res[!is_zero] <- rlnorm(sum(!is_zero), meanlog = meanlog, sdlog = sdlog)

  return(res)
}

#' @rdname zilnorm
#' @importFrom RTMB pnorm
#' @export
plnorm <- function(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE) {
  if(!ad_context()) {
    # ensure sdlog > 0
    if (any(sdlog <= 0)) stop("sdlog must be > 0")
  }

  eps <- .Machine$double.xmin
  ipos <- ispos_strict(q)        # 1 for q > 0, 0 otherwise
  q_safe <- q * ipos + eps       # = q (for q>0, +eps) or eps (for q<=0)
  z <- (log(q_safe) - meanlog) / sdlog

  p <- ipos * pnorm(z)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}

