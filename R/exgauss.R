#' Exponentially modified Gaussian distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the exponentially modified Gaussian distribution.
#'
#' @details
#' This implementation of \code{dexgauss} and \code{pexgauss} allows for automatic differentiation with \code{RTMB}.
#' \code{qexgauss} inverts \code{pexgauss} numerically, as the exponentially modified
#' Gaussian has no closed-form quantile function.
#' The parameterisation follows the \code{exGAUS} family of the \code{gamlss.dist}
#' package, with the exponential rate \eqn{\lambda = 1/\nu}.
#'
#' If \eqn{X \sim N(\mu, \sigma^2)} and \eqn{Y \sim \text{Exp}(\lambda)}, then
#' \eqn{Z = X + Y} follows the exponentially modified Gaussian distribution with parameters \eqn{\mu}, \eqn{\sigma}, and \eqn{\lambda}.
#'
#' The density is
#' \deqn{f(x;\,\mu,\sigma,\lambda) = \lambda \exp\!\Bigl(\lambda\mu + \tfrac{\lambda^2\sigma^2}{2} - \lambda x\Bigr)\, \Phi\!\left(\frac{x - \mu - \lambda\sigma^2}{\sigma}\right),}
#' where \eqn{\Phi} is the standard normal CDF.
#'
#' @references
#' Cousineau, D., Brown, S. and Heathcote, A. (2004) Fitting distributions using maximum likelihood: Methods and packages.
#' Behavior Research Methods, Instruments, & Computers, 36, 742-756.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [jsu], [skewnorm]
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu mean parameter of the Gaussian part
#' @param sigma standard deviation parameter of the Gaussian part, must be positive.
#' @param lambda rate parameter of the exponential part, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dexgauss} gives the density, \code{pexgauss} gives the distribution function, \code{qexgauss} gives the quantile function, and \code{rexgauss} generates random deviates.
#'
#' @examples
#' x <- rexgauss(1, 1, 2, 2)
#' d <- dexgauss(x, 1, 2, 2)
#' p <- pexgauss(x, 1, 2, 2)
#' q <- qexgauss(p, 1, 2, 2)
#' @name exgauss
NULL

#' @rdname exgauss
#' @export
#' @importFrom RTMB dnorm pnorm
dexgauss <- function(x, mu = 0, sigma = 1, lambda = 1, log = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/exGAUS.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure sigma > 0, lambda > 0
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dexgauss", x=x, mu=mu, sigma=sigma, lambda=lambda, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dexgauss", x=x, mu=mu, sigma=sigma, lambda=lambda, log=log))
  }

  nu <- 1 / lambda

  z <- x - mu - ((sigma * sigma) * lambda)

  nu_gr <- greater(nu, 0.05 * sigma) # numerical stability

  logdens <- nu_gr * as.finite.neg(- log(nu) - (z + ((sigma * sigma) / (2 * nu))) / nu + log(1e-300 + RTMB::pnorm(z / sigma))) +
    (1-nu_gr) * RTMB::dnorm(x, mean = mu, sd = sigma, log = TRUE)

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname exgauss
#' @export
#' @importFrom RTMB pnorm
pexgauss <- function(q, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/exGAUS.R
  # and modified to allow for automatic differentiation

  if (!ad_context()) {
    # ensure sigma > 0, lambda > 0
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda must be > 0")
  }

  nu <- 1 / lambda

  z <- q - mu - (sigma^2 / nu)

  nu_gr <- greater(nu, 0.05 * sigma) # numerical stability

  p <- nu_gr * RTMB::pnorm((q-mu)/sigma) - RTMB::pnorm(z/sigma) * exp(((mu+(sigma^2/nu))^2-(mu^2)-2*q*((sigma^2)/nu))/(2*sigma^2)) +
    (1 - nu_gr) * RTMB::pnorm(q, mean = mu, sd = sigma)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname exgauss
#' @export
qexgauss <- function(p, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/exGAUS.R
  # there is no closed form, so the cdf is bracketed outwards from mu in steps
  # of sigma and then inverted with uniroot

  if (!ad_context()) {
    # ensure sigma > 0, lambda > 0
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda must be > 0")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  lp <- max(lengths(list(p, mu, sigma, lambda)))
  p <- rep_len(p, lp); mu <- rep_len(mu, lp)
  sigma <- rep_len(sigma, lp); lambda <- rep_len(lambda, lp)

  q <- numeric(lp)
  for (i in seq_len(lp)) {

    if (p[i] <= 0) { q[i] <- -Inf; next }
    if (p[i] >= 1) { q[i] <-  Inf; next }

    h <- function(z) pexgauss(z, mu = mu[i], sigma = sigma[i], lambda = lambda[i])

    # step outwards from mu in multiples of sigma until the root is bracketed;
    # s = +1 searches upwards, s = -1 downwards
    s <- if (h(mu[i]) < p[i]) 1 else -1
    k <- 1
    repeat {
      edge <- mu[i] + s * k * sigma[i]
      if (s * (h(edge) - p[i]) >= 0 || k >= 1e6) break
      k <- k + 1
    }

    q[i] <- stats::uniroot(function(z) h(z) - p[i], sort(c(mu[i], edge)),
                           tol = .Machine$double.eps^0.75)$root
  }

  return(q)
}
#' @rdname exgauss
#' @export
#' @importFrom stats rnorm rexp
rexgauss <- function(n, mu = 0, sigma = 1, lambda = 1) {
  # ensure sigma > 0, lambda > 0
  if (any(sigma <= 0)) stop("sigma must be > 0")
  if (any(lambda <= 0)) stop("lambda must be > 0")

  # Generate n random values from the exponentially modified Gaussian distribution
  rnorm(n, mu, sigma) + rexp(n, rate = lambda)
}
