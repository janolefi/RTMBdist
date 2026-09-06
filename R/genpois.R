#' Generalised Poisson distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the generalised Poisson distribution.
#'
#' @details
#' This implementation of \code{dgenpois} allows for automatic differentiation with \code{RTMB}.
#' The parameterisation follows the \code{GPO} family of the \code{gamlss.dist} package.
#'
#' The distribution has mean \eqn{\lambda} and variance \eqn{\lambda(1 + \phi \lambda)^2}.
#' For \eqn{\phi = 0} it reduces to the Poisson distribution, however \eqn{\phi} must be strictly positive here.
#'
#' \deqn{P(X = x;\,\lambda,\phi) = \frac{\lambda\,(1+\phi x)^{x-1}\,e^{-\lambda(1+\phi x)/(1+\phi\lambda)}}{(1+\phi\lambda)^x\, x!}, \quad x = 0, 1, 2, \ldots}
#'
#' @references
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [nbinom2], [zipois], [bell]
#'
#' @param x,q integer vector of counts
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param lambda vector of positive means
#' @param phi vector of non-negative dispersion parameters
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' @param max.value a constant, set to the default value of 10000 for how far the algorithm should look for \code{q}.
#'
#' @return
#' \code{dgenpois} gives the probability mass function, \code{pgenpois} gives the distribution function, \code{qgenpois} gives the quantile function, and \code{rgenpois} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rgenpois(1, 2, 3)
#' d <- dgenpois(x, 2, 3)
#' p <- pgenpois(x, 2, 3)
#' q <- qgenpois(p, 2, 3)
#' @name genpois
NULL
#' @rdname genpois
#' @export
#' @import RTMB
dgenpois <- function(x, lambda = 1, phi = 1, log = FALSE) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/GPO.R
  # and modified to allow for automatic differentiation

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure lambda, phi > 0
    if (any(lambda <= 0)) stop("lambda must be > 0")
    if (any(phi <= 0)) stop("phi must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dgenpois", x = x, lambda = lambda, phi = phi, log=log))
  }
  if(inherits(x, "osa")) {
    # OSA works via method = "oneStepGeneric" (density only); method = "cdf"
    # is unavailable because pgenpois sums the pmf in an R loop and cannot be taped
    return(dGenericOSA("dgenpois", x = x, lambda = lambda, phi = phi, log=log))
  }

  phi_lambda <- phi * lambda
  phi_x <- phi * x

  logdens <- x * (log(lambda) - log1p(phi_lambda)) +
    (x-1) * log1p(phi_x) - lgamma(x + 1) -
    (lambda * (phi_x + 1)) / (phi_lambda + 1)

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname genpois
#' @export
pgenpois <- function(q, lambda = 1, phi = 1, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure lambda, phi > 0
    if (any(lambda <= 0)) stop("lambda must be > 0")
    if (any(phi <= 0)) stop("phi must be > 0")
    # Check q is integer >= 0
    if(any(q < 0 | q != floor(q))) {
      stop("q must be a non-negative integer")
    }
  }

  # summing the pmf over 0:q is exactly what gamlss.dist::pGPO does
  pgenpois.ad(q=q, lambda=lambda, phi=phi, lower.tail=lower.tail, log.p=log.p)
}
pgenpois.ad <- function(q, lambda = 1, phi = 1, lower.tail = TRUE, log.p = FALSE){

  # summing the pmf over 0:q is the same approach taken by
  # https://github.com/gamlss-dev/gamlss.dist/blob/main/R/GPO.R (pGPO),
  # written here so that it also works under automatic differentiation

  # a single q with vectorised parameters must still give one value per parameter
  n <- max(length(q), length(lambda), length(phi))
  if(length(q) != n) q <- rep(q, length.out = n)
  if(length(lambda) != n) lambda <- rep(lambda, length.out = n)
  if(length(phi) != n) phi <- rep(phi, length.out = n)

  p <- rep(0, n)

  for(i in seq_len(n)) {
    x <- 0:q[i]
    p[i] <- sum(dgenpois(x, lambda[i], phi[i]))
  }

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)

  return(p)
}
#' @rdname genpois
#' @export
#' @usage qgenpois(p, lambda = 1, phi = 1,
#'          lower.tail = TRUE, log.p = FALSE, max.value = 10000)
qgenpois <- function(p, lambda = 1, phi = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1e4) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/GPO.R

  if(!ad_context()) {
    # ensure lambda, phi > 0
    if (any(lambda <= 0)) stop("lambda must be > 0")
    if (any(phi <= 0)) stop("phi must be > 0")
  }

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p

  if(!ad_context()) {
    # Check p is in [0,1]
    if(any(p < 0 | p > 1)) {
      stop("p must be in [0,1]")
    }
  }

  # a single p with vectorised parameters must still give one value per parameter
  ly <- max(lengths(list(p, lambda, phi)))
  p <- rep_len(p, ly)
  lambda <- rep_len(lambda, ly)
  phi <- rep_len(phi, ly)

  support <- 0:max.value
  q <- numeric(ly)

  # the cdf only has to be rebuilt when the parameters change; `built` is the
  # index the cached cdf belongs to (0 = nothing cached yet)
  cdf <- NULL
  built <- 0L

  for(i in seq_len(ly)) {
    if(p[i] + 1e-09 >= 1) {
      q[i] <- Inf
      next
    }
    if(built == 0L || lambda[i] != lambda[built] || phi[i] != phi[built]) {
      cdf <- cumsum(dgenpois(support, lambda[i], phi[i]))
      built <- i
    }
    j <- match(TRUE, p[i] <= cdf)
    # gamlss.dist returns max.value when the search does not reach p
    q[i] <- if(is.na(j)) max.value else support[j]
  }

  return(q)
}
#' @rdname genpois
#' @export
#' @importFrom stats runif
rgenpois <- function(n, lambda = 1, phi = 1, max.value = 1e4) {

  # taken from https://github.com/gamlss-dev/gamlss.dist/blob/main/R/GPO.R

  if(!ad_context()) {
    # ensure lambda, phi > 0
    if (any(lambda <= 0)) stop("lambda must be > 0")
    if (any(phi <= 0)) stop("phi must be > 0")
  }

  n <- ceiling(n)
  p <- runif(n)

  qgenpois(p, lambda = lambda, phi = phi, max.value = max.value)
}

