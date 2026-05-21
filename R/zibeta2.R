#' Reparameterised zero-inflated beta distribution
#'
#' Density, distribution function, and random generation for
#' the zero-inflated beta distribution reparameterised in terms of mean and concentration.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' Uses the same density as \code{zibeta} with \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}:
#' \deqn{f(x;\,\mu,\phi,p_0) = p_0\,\mathbf{1}[x=0] + (1-p_0)\,f_{\mathrm{Beta}}(x;\,\mu\phi,\,(1-\mu)\phi)\,\mathbf{1}[x\in(0,1)].}
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mu mean parameter, must be in the interval from 0 to 1.
#' @param phi concentration parameter, must be positive.
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dzibeta2} gives the density, \code{pzibeta2} gives the distribution function, and \code{rzibeta2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzibeta2(1, 0.5, 1, 0.5)
#' d <- dzibeta2(x, 0.5, 1, 0.5)
#' p <- pzibeta2(x, 0.5, 1, 0.5)
#' @name zibeta2
NULL
#' @rdname zibeta2
#' @export
dzibeta2 <- function(x, mu, phi, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
    # ensure zeroprob in [0,1]
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dzibeta2", x=x, mu=mu, phi=phi, zeroprob=zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzibeta2", x=x, mu=mu, phi=phi, zeroprob=zeroprob, log=log))
  }

  # parameter transformation
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  dzibeta(x, shape1 = shape1, shape2 = shape2, zeroprob = zeroprob, log = log)
}
#' @rdname zibeta2
#' @export
pzibeta2 <- function(q, mu, phi, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
    # ensure zeroprob in [0,1]
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # parameter transformation
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  pzibeta(q, shape1 = shape1, shape2 = shape2, zeroprob = zeroprob,
          lower.tail = lower.tail, log.p = log.p)
}
#' @rdname zibeta2
#' @export
rzibeta2 <- function(n, mu, phi, zeroprob = 0) {

  if(!ad_context()) {
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
    # ensure zeroprob in [0,1]
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # parameter transformation
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  rzibeta(n, shape1 = shape1, shape2 = shape2, zeroprob = zeroprob)
}


