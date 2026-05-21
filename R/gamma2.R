#' Reparameterised gamma distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the gamma distribution reparameterised in terms of mean and standard deviation.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' Reparameterises the gamma distribution via \eqn{\text{shape} = \mu^2/s^2} and \eqn{\text{scale} = s^2/\mu}:
#' \deqn{f(x;\,\mu, s) = \frac{x^{\mu^2/s^2 - 1} \exp(-x s^2/\mu^2)}{(s^2/\mu)^{\mu^2/s^2}\,\Gamma(\mu^2/s^2)}, \quad x > 0.}
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mean mean parameter, must be positive.
#' @param sd standard deviation parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dgamma2} gives the density, \code{pgamma2} gives the distribution function, \code{qgamma2} gives the quantile function, and \code{rgamma2} generates random deviates.
#'
#' @examples
#' x <- rgamma2(1)
#' d <- dgamma2(x)
#' p <- pgamma2(x)
#' q <- qgamma2(p)
#' @name gamma2
NULL

#' @rdname gamma2
#' @export
#' @importFrom RTMB dgamma
dgamma2 = function(x, mean = 1, sd = 1, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure mean, sd > 0
    if (any(mean <= 0)) stop("mean must be strictly positive.")
    if (any(sd <= 0)) stop("sd must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dgamma2", x=x, mean = mean, sd = sd, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dgamma2", x=x, mean = mean, sd = sd, log=log))
  }

  # parameter transformation
  shape <- mean^2 / sd^2
  scale <- sd^2 / mean
  dgamma(x = x, shape = shape, scale = scale, log = log)
}

#' @rdname gamma2
#' @export
#' @importFrom RTMB pgamma
pgamma2 = function(q, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure mean, sd > 0
    if (any(mean <= 0)) stop("mean must be strictly positive.")
    if (any(sd <= 0)) stop("sd must be strictly positive.")
  }

  # parameter transformation
  shape <- mean^2 / sd^2
  scale <- sd^2 / mean
  p <- pgamma(q = q, shape = shape, scale = scale)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}

#' @rdname gamma2
#' @export
#' @importFrom RTMB qgamma
qgamma2 = function(p, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure mean, sd > 0
    if (any(mean <= 0)) stop("mean must be strictly positive.")
    if (any(sd <= 0)) stop("sd must be strictly positive.")
  }

  # parameter transformation
  shape <- mean^2 / sd^2
  scale <- sd^2 / mean
  qgamma(p = p, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gamma2
#' @export
#' @importFrom stats rgamma
rgamma2 = function(n, mean = 1, sd = 1) {
  if(!ad_context()) {
    # ensure mean, sd > 0
    if (any(mean <= 0)) stop("mean must be strictly positive.")
    if (any(sd <= 0)) stop("sd must be strictly positive.")
  }

  # parameter transformation
  shape = mean^2 / sd^2
  scale = sd^2 / mean
  rgamma(n = n, shape = shape, scale = scale)
}
