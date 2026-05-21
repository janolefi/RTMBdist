#' Student t distribution with location and scale
#'
#' Density, distribution function, quantile function, and random generation for
#' the t distribution with location and scale parameters.
#'
#' @details
#' This implementation of \code{dt2} allows for automatic differentiation with \code{RTMB}.
#'
#' \deqn{f(x;\,\mu,\sigma,\nu) = \frac{1}{\sigma}\,f_t\!\left(\frac{x-\mu}{\sigma};\,\nu\right),}
#' where \eqn{f_t(\cdot;\nu)} is the Student-\eqn{t} PDF with \eqn{\nu} degrees of freedom.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mu location parameter
#' @param sigma scale parameter, must be positive.
#' @param df degrees of freedom, must be positive.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{dt2} gives the density, \code{pt2} gives the distribution function, \code{qt2} gives the quantile function, and \code{rt2} generates random deviates.
#'
#' @examples
#' x <- rt2(1, 1, 2, 5)
#' d <- dt2(x, 1, 2, 5)
#' p <- pt2(x, 1, 2, 5)
#' q <- qt2(p, 1, 2, 5)
#' @name t2
NULL
#' @rdname t2
#' @export
#' @importFrom RTMB dt
dt2 = function(x, mu, sigma, df, log = FALSE){

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure sigma > 0
    if (any(sigma <= 0)) stop("sigma must be strictly positive.")
    # ensure df > 0
    if (any(df <= 0)) stop("df must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dt2", x=x, mu=mu, sigma=sigma, df=df, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dt2", x=x, mu=mu, sigma=sigma, df=df, log=log))
  }

  z <- (x - mu) / sigma
  logdens <- RTMB::dt(z, df, log = TRUE) - log(sigma)

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname t2
#' @export
#' @importFrom stats pt
pt2 <- function(q, mu, sigma, df, lower.tail = TRUE, log.p = FALSE){
  z <- (q - mu) / sigma
  p <- pt(z, df)
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}
#' @rdname t2
#' @export
#' @importFrom stats rt
rt2 <- function(n, mu, sigma, df) {
  z <- stats::rt(n, df)
  return(mu + sigma * z)
}
#' @rdname t2
#' @export
#' @importFrom stats qt
qt2 <- function(p, mu, sigma, df, lower.tail = TRUE, log.p = FALSE){
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  z <- stats::qt(p, df)
  return(mu + sigma * z)
}

#' @rdname t2
#' @export
#' @importFrom RTMB pbeta
pt <- function(q, df) {
  # AD compatible version of pt - slightly sketchy but works

  x <- df / (df + q^2)
  q <- q + 1e-8 # avoid numerical issues with q = 0
  # if skew is exactly zero, q will be zero and the gradient will be exactly zero

  val <- RTMB::pbeta(x, df / 2, 0.5) / 2
  test <- 0.5 * (sign(q) + 1)  # test if q > 0
  test * (1 - val) + (1 - test) * val
}
