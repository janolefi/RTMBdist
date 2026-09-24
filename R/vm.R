#' von Mises distribution
#'
#' Density, distribution function, and random generation for the von Mises distribution.
#'
#' @details
#' This implementation of \code{dvm} allows for automatic differentiation with \code{RTMB}.
#' \code{rvm} and \code{pvm} are simply wrappers of the corresponding functions from \code{circular}.
#' When called during AD taping, \code{pvm} instead integrates the density numerically with the AD-compatible \code{\link[RTMB:ADintegrate]{integrate}} of \code{RTMB}, which makes it AD-compatible in \code{q}, \code{mu}, \code{kappa} and \code{from}. The \code{tol} argument is then ignored.
#'
#' \deqn{f(x;\,\mu,\kappa) = \frac{\exp(\kappa\cos(x-\mu))}{2\pi\, I_0(\kappa)},}
#' where \eqn{I_0} is the modified Bessel function of the first kind of order 0.
#'
#' A circular distribution has no smallest angle, so its distribution function
#' depends on where the circle is cut open. By default, \code{pvm} cuts it at the
#' antipode of the mean direction, \eqn{\mu - \pi}, so that \eqn{F(\mu) = 1/2}.
#' A different origin is set with \code{from}. A fixed \code{from} is needed
#' whenever the distribution function is averaged over different values of
#' \code{mu}, for example over the states of a hidden Markov model or over a
#' random effect: with the default, each value of \code{mu} cuts the circle at a
#' different place, and the average of these distribution functions is not the
#' distribution function of the mixture.
#'
#' \strong{OSA residuals:} \code{dvm} supports one-step-ahead (OSA) quantile
#' residuals via \code{RTMB::\link[RTMB]{oneStepPredict}}. For the methods based
#' on the distribution function, such as \code{method = "cdf"}, the circle is cut
#' at the fixed origin \eqn{-\pi}, i.e. the residuals are based on
#' \code{pvm(x, mu, kappa, from = -pi)} rather than the default origin
#' \eqn{\mu - \pi}. OSA residuals are computed from the predictive distribution
#' function, which averages the distribution function over hidden states or
#' random effects, and this is only valid with an origin that does not depend on
#' \code{mu} (see above). Hence the residuals are valid for all models, but their
#' interpretation depends on the data: for turning angles, with \code{mu} close
#' to 0, the cut at \eqn{\pm\pi} corresponds to a reversal and the residuals
#' increase with the turning angle. For directions with \code{mu} far from 0,
#' angles close to \eqn{\pm\pi} can give large residuals of either sign, even
#' when they are close to the mean direction. For
#' \code{method = "oneStepGeneric"}, set \code{range = c(-pi, pi)} in
#' \code{oneStepPredict()}, so that the density is integrated from the same origin.
#'
#' @param x,q vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param kappa non-negative numeric value for the concentration parameter of the distribution.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of random values to return.
#' @param tol the precision in evaluating the distribution function, ignored in AD context.
#' @param from value from which the integration for CDF starts. If \code{NULL}, is set to \code{mu - pi}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param log.p logical; if \code{TRUE}, probabilities are returned as \eqn{\log(p)}.
#' @param wrap logical; if \code{TRUE}, generated angles are wrapped to the interval from -pi to pi.
#'
#' @return \code{dvm} gives the density, \code{pvm} gives the distribution function, and \code{rvm} generates random deviates.
#'
#' @seealso [wrpcauchy]; [cjw()] and [cfold()] for circular-linear copulas joining turning angles and step lengths.
#'
#' @examples
#' set.seed(1)
#' x <- rvm(10, 0, 1)
#' d <- dvm(x, 0, 1)
#' p <- pvm(x, 0, 1)
#' @name vm
NULL

#' @rdname vm
#' @export
#' @importFrom RTMB besselI
dvm = function(x, mu = 0, kappa = 1, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure kappa >= 0
    if (any(kappa < 0)) stop("kappa must be non-negative.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dvm", x = x, mu = mu, kappa = kappa, log=log))
  }
  if(inherits(x, "osa")) {
    # the circle is cut at the fixed origin -pi, see the OSA section of the documentation
    return(dGenericOSA("dvm_osa", x = x, mu = mu, kappa = kappa, log = log))
  }

  # stable calculation of log(besselI(kappa, 0))
  logI0 <- log(RTMB::besselI(kappa, 0, expon.scaled = TRUE)) + kappa

  logdens <- -log(2 * pi) - logI0 + kappa * cos(x - mu)

  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
  }
}

# density and distribution function behind the OSA residuals of dvm, with the circle
# cut at -pi, found by RTMB's dGenericOSA via their names
dvm_osa <- function(x, mu, kappa, log = FALSE) dvm(x, mu, kappa, log = log)
pvm_osa <- function(q, mu, kappa) pvm(q, mu, kappa, from = -pi)

#' @rdname vm
#' @export
#' @importFrom circular pvonmises
pvm = function(q, mu = 0, kappa = 1, from = NULL, tol = 1e-20,
               lower.tail = TRUE, log.p = FALSE) {

  if (ad_context()) {
    # circular::pvonmises is not AD-compatible, hence integrate the density numerically.
    # The angle from mu is wrapped to [-pi, pi) with floor(), which is re-evaluated with the
    # tape (%% is not for advectors), so mu may be a parameter and angles may cross mu +- pi.
    wrap <- function(a) a - 2 * pi * floor((a + pi) / (2 * pi))
    cdf <- function(a, ...) numerical_cdf(dvm, wrap(a - mu), list(mu = 0, kappa = kappa),
                                         centre = 0, lower = -pi, upper = pi, ...)
    if (is.null(from)) return(cdf(q, lower.tail = lower.tail, log.p = log.p))
    # fixed origin: (F(q) - F(from)) mod 1, the wrap again written with floor()
    probs <- cdf(q) - cdf(from)
    probs <- probs - floor(probs)
    if (!lower.tail) probs <- 1 - probs
    if (log.p) probs <- log(probs)
    return(probs)
  }
  # NA handling
  ind = which(!is.na(q))

  if(is.matrix(mu)){
    mu = mu[ind,]
  }
  if(is.matrix(kappa)){
    kappa = kappa[ind,]
  }

  probs = numeric(length(q))

  suppressWarnings(
    probs[ind] <- pvonmises(q[ind], mu, kappa, from = from, tol = tol)
  )

  probs[-ind] = NA

  probs <- as.numeric(probs)
  if (!lower.tail) probs <- 1 - probs
  if (log.p) probs <- log(probs)
  probs
}

#' @rdname vm
#' @export
#' @importFrom circular rvonmises
rvm = function(n, mu = 0, kappa = 1, wrap = TRUE) {
  suppressWarnings(
    angles <- as.numeric(rvonmises(n, mu, kappa))
  )

  # if generated angels should be wrapped, i.e. mapped to interval [-pi, pi], do so
  if(wrap){
    angles = (angles + pi) %% (2 * pi) - pi
  }
  angles
}
