#' Skewed students t distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the skew t distribution (type 2).
#'
#' @details
#' This corresponds to the skew t type 2 distribution in GAMLSS (\code{\link[gamlss.dist]{ST2}}), see pp. 411-412 of Rigby et al. (2019) and the version implemented in the \code{sn} package.
#' This implementation of \code{dskewt} allows for automatic differentiation with \code{RTMB} while the other functions are imported from the \code{sn} package.
#' See \code{sn::\link[sn]{dst}} for more details.
#'
#' \strong{Caution:} In a numerial optimisation, the \code{skew} parameter should NEVER be initialised with exactly zero.
#' This will cause the initial and all subsequent derivatives to be exactly zero and hence the parameter will remain at its initial value.
#'
#' @seealso [skewt2], [skewnorm], [skewnorm2]
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mu location parameter
#' @param sigma scale parameter, must be positive.
#' @param skew skewness parameter, can be positive or negative.
#' @param df degrees of freedom, must be positive.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param tol a scalar value which regulates the accuracy of the result of qsn, measured on the probability scale.
#' @param method an integer value between 0 and 5 which selects the computing method; see ‘Details’ in the \code{\link[sn]{pst}} documentation below for the meaning of these values. If method=0 (default value), an automatic choice is made among the four actual computing methods, depending on the other arguments.
#'
#' @return
#' \code{dskewt} gives the density, \code{pskewt} gives the distribution function, \code{qskewt} gives the quantile function, and \code{rskewt} generates random deviates.
#'
#' @examples
#' x <- rskewt(1, 1, 2, 5, 2)
#' d <- dskewt(x, 1, 2, 5, 2)
#' p <- pskewt(x, 1, 2, 5, 2)
#' q <- qskewt(p, 1, 2, 5, 2)
#' @name skewt
NULL
#' @rdname skewt
#' @export
#' @importFrom RTMB dt
dskewt <- function(x, mu = 0, sigma = 1, skew = 0, df = 1e2, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure sigma, df > 0
    if (sigma <= 0) stop("sigma must be strictly positive.")
    if (df <= 0) stop("df must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dskewt", x=x, mu=mu, sigma=sigma, skew=skew, df=df, log=log))
  }
  if(inherits(x, "osa")) {
    # return(dGenericOSA("dskewt", x=x, mu=mu, sigma=sigma, skew=skew, df=df, log=log))
    stop("Currently, skew t does not support OSA residuals.")
  }

  z <- (x - mu) / sigma
  lambda <- (df + 1) / (df + z * z)
  omega <- skew * sqrt(lambda) * z

  pdf <- RTMB::dt(z, df, log = TRUE)
  cdf <- log(1e-300 + pt(omega, df + 1))

  logdens <- log(2) - log(sigma) + pdf + cdf

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname skewt
#' @export
#' @usage
#' pskewt(q, mu = 0, sigma = 1, skew = 0, df = 100,
#'        method = 0, lower.tail = TRUE, log.p = FALSE)
#' @importFrom sn pst
pskewt <- function(q, mu = 0, sigma = 1, skew = 0, df = 1e2, method = 0, lower.tail = TRUE, log.p = FALSE) {
  # ensure sigma, df > 0
  # if (sigma <= 0) stop("sigma must be strictly positive.")
  # if (df <= 0) stop("df must be strictly positive.")

  pst(q, xi=mu, omega=sigma, alpha=skew, nu=df, method=method, lower.tail=lower.tail, log.p=log.p)
}

#' @rdname skewt
#' @export
#' @usage
#' qskewt(p, mu = 0, sigma = 1, skew = 0, df = 100,
#'        tol = 1e-8, method = 0)
#' @importFrom sn qst
qskewt <- function(p, mu = 0, sigma = 1, skew = 0, df = 1e2, tol = 1e-8, method = 0) {
  # ensure sigma, df > 0
  if (sigma <= 0) stop("sigma must be strictly positive.")
  if (df <= 0) stop("df must be strictly positive.")
  qst(p, xi=mu, omega=sigma, alpha=skew, nu=df, tol=tol, method=method)
}

#' @rdname skewt
#' @export
#' @importFrom sn rst
rskewt <- function(n, mu = 0, sigma = 1, skew = 0, df = 1e2) {
  # ensure sigma, df > 0
  if (sigma <= 0) stop("sigma must be strictly positive.")
  if (df <= 0) stop("df must be strictly positive.")

  as.numeric(rst(n, xi=mu, omega=sigma, alpha=skew, nu=df))
}
