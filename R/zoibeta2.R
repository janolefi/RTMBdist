#' Reparameterised zero- and one-inflated beta distribution
#'
#' Density, distribution function, and random generation for
#' the zero-one-inflated beta distribution reparameterised in terms of mean and concentration.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return.
#' @param mu mean parameter, must be in the interval from 0 to 1.
#' @param phi concentration parameter, must be positive.
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param oneprob one-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzoibeta2} gives the density, \code{pzoibeta2} gives the distribution function, and \code{rzoibeta2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzoibeta2(1, 0.6, 2, 0.2, 0.3)
#' d <- dzoibeta2(x, 0.6, 2, 0.2, 0.3)
#' p <- pzoibeta2(x, 0.6, 2, 0.2, 0.3)
#' @name zoibeta2
NULL
#' @rdname zoibeta2
#' @export
dzoibeta2 <- function(x, mu, phi, zeroprob = 0, oneprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dzoibeta2", x=x, mu=mu, phi=phi,
                       zeroprob=zeroprob, oneprob=oneprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzoibeta2", x=x, mu=mu, phi=phi,
                       zeroprob=zeroprob, oneprob=oneprob, log=log))
  }

  # parameter transformation
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  dzoibeta(x, shape1 = shape1, shape2 = shape2,
           zeroprob = zeroprob, oneprob = oneprob, log = log)
}
#' @rdname zoibeta2
#' @usage pzoibeta2(q, mu, phi, zeroprob = 0, oneprob = 0,
#'          lower.tail = TRUE, log.p = FALSE)
#' @export
pzoibeta2 <- function(q, mu, phi, zeroprob = 0, oneprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure mu in [0,1]
    if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
    # ensure phi > 0
    if (any(phi <= 0)) stop("phi must be strictly positive.")
  }

  # parameter transformation
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  pzoibeta(q, shape1 = shape1, shape2 = shape2, zeroprob = zeroprob, oneprob = oneprob,
          lower.tail = lower.tail, log.p = log.p)
}
#' @rdname zoibeta2
#' @export
rzoibeta2 <- function(n, mu, phi, zeroprob = 0, oneprob = 0) {

  # ensure mu in [0,1]
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
  # ensure phi > 0
  if (any(phi <= 0)) stop("phi must be strictly positive.")

  # parameter transformation
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  rzoibeta(n, shape1 = shape1, shape2 = shape2,
           zeroprob = zeroprob, oneprob = oneprob)
}
