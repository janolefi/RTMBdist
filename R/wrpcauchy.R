#' wrapped Cauchy distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the wrapped Cauchy distribution.
#'
#' @details
#' \code{dwrpcauchy} and \code{pwrpcauchy} allow for automatic differentiation with \code{RTMB}.
#'
#' \deqn{f(x;\,\mu,\rho) = \frac{1}{2\pi}\cdot\frac{1-\rho^2}{1 + \rho^2 - 2\rho\cos(x-\mu)}.}
#'
#' A circular distribution has no smallest angle, so its distribution function
#' depends on where the circle is cut open. \code{pwrpcauchy} cuts it at the
#' antipode of the mean direction, \eqn{\mu - \pi}, and is then available in
#' closed form:
#' \deqn{F(q) = P(\mu - \pi < X \le q) = \frac{1}{2} + \frac{1}{\pi}\arctan\left(\frac{1+\rho}{1-\rho}\tan\frac{q-\mu}{2}\right), \quad \mu - \pi < q \le \mu + \pi,}
#' so that \eqn{F(\mu) = 1/2}. This is the default, \code{from = NULL}, and
#' the same origin as the default of \code{\link{pvm}}.
#'
#' Angles outside \eqn{(\mu - \pi, \mu + \pi]} are wrapped onto this interval
#' first, so \code{pwrpcauchy} is \eqn{2\pi}-periodic in \code{q}. As a
#' consequence, the default origin moves with \code{mu}: for angles on
#' \eqn{[-\pi, \pi]} and \eqn{\mu \neq 0}, \code{pwrpcauchy} is not monotone
#' over that range but drops from 1 back to 0 at \eqn{\mu \pm \pi}.
#'
#' A different origin is set with \code{from}, giving
#' \eqn{P(\mathrm{from} < X \le q) = (F(q) - F(\mathrm{from})) \bmod 1} for
#' \eqn{\mathrm{from} < q \le \mathrm{from} + 2\pi}. A fixed \code{from} is
#' needed whenever the distribution function is averaged over different values
#' of \code{mu}, for example over the states of a hidden Markov model or over a
#' random effect: with the default, each value of \code{mu} cuts the circle at a
#' different place, and the average of these distribution functions is not the
#' distribution function of the mixture.
#'
#' \code{qwrpcauchy} is the inverse of \code{pwrpcauchy} and returns angles in
#' \eqn{[\mathrm{from}, \mathrm{from} + 2\pi]}, by default
#' \eqn{[\mu - \pi, \mu + \pi]}, not wrapped to \eqn{[-\pi, \pi]}. The latter
#' also holds for \code{rwrpcauchy} with \code{wrap = FALSE}.
#'
#' @seealso [vm]; [cjw()] and [cfold()] for circular-linear copulas joining turning angles and step lengths.
#'
#' @param x,q vector of angles measured in radians at which to evaluate the density or distribution function.
#' @param p vector of probabilities.
#' @param mu mean direction of the distribution measured in radians.
#' @param rho concentration parameter of the distribution, must be in the interval from 0 to 1.
#' @param from origin, in radians, at which the circle is cut open for the distribution and quantile functions. If \code{NULL} (default), it is set to \code{mu - pi}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param n number of random values to return.
#' @param wrap logical; if \code{TRUE}, generated angles are wrapped to the interval from -pi to pi, otherwise they lie in the interval from \code{mu - pi} to \code{mu + pi}.
#'
#' @return \code{dwrpcauchy} gives the density, \code{pwrpcauchy} gives the distribution function, \code{qwrpcauchy} gives the quantile function, and \code{rwrpcauchy} generates random deviates.
#'
#' @examples
#' set.seed(1)
#' x <- rwrpcauchy(10, 0, 0.5)
#' d <- dwrpcauchy(x, 0, 0.5)
#' p <- pwrpcauchy(x, 0, 0.5)
#' q <- qwrpcauchy(p, 0, 0.5)
#' @name wrpcauchy
NULL

#' @rdname wrpcauchy
#' @export
dwrpcauchy <- function(x, mu = 0, rho, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure rho in [0, 1)
    if (any(rho < 0) || any(rho >= 1)) stop("rho must be in the interval [0, 1).")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dwrpcauchy", x = x, mu = mu, rho = rho, log=log))
  }
  if(inherits(x, "osa")) {
    stop("Wrapped cauchy does not support OSA residuals.")
  }
  rho_sq <- rho * rho

  logdens <- - log(2 * pi) +
    log1p(-rho_sq) -
    log1p(rho_sq - 2 * rho * cos(x - mu))

  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
  }
}

#' @rdname wrpcauchy
#' @export
pwrpcauchy <- function(q, mu = 0, rho, from = NULL, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(rho < 0) || any(rho >= 1)) stop("rho must be in the interval [0, 1).")
  }

  # tan() has period pi in (q - mu) / 2, which wraps q onto (mu - pi, mu + pi].
  # At the cut itself, tan(pi/2) is merely large in double precision, so the
  # value stays finite and so does its derivative.
  z <- atan((1 + rho) / (1 - rho) * tan((q - mu) / 2)) / pi

  if (is.null(from)) {
    # the upper tail is computed directly rather than as 1 - p
    if (lower.tail) {
      p <- 0.5 + z
    } else {
      p <- 0.5 - z
    }
  } else {
    # probability from `from` to q counterclockwise; %% has derivative 1
    zf <- atan((1 + rho) / (1 - rho) * tan((from - mu) / 2)) / pi
    p <- (z - zf) %% 1
    if (!lower.tail) p <- 1 - p
  }
  if (log.p) p <- log(p)
  return(p)
}

#' @rdname wrpcauchy
#' @export
qwrpcauchy <- function(p, mu = 0, rho, from = NULL, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(rho < 0) || any(rho >= 1)) stop("rho must be in the interval [0, 1).")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  # quantile function for the origin mu - pi, returning angles in [mu - pi, mu + pi];
  # at p = 0 and p = 1, tan(-+pi/2) is large enough for atan to return -+pi/2,
  # so the end points need no special treatment
  q0 <- function(p) mu + 2 * atan((1 - rho) / (1 + rho) * tan(pi * (p - 0.5)))

  if (is.null(from)) return(q0(p))

  # g is the probability from mu - pi, running past 1 once the quantile passes
  # mu + pi, in which case it continues one turn further
  pf <- pwrpcauchy(from, mu, rho)
  g <- p + pf
  w <- greater(g, 1)
  q <- q0(g - w) + 2 * pi * w
  # q starts at q0(pf), which is from wrapped onto [mu - pi, mu + pi]; shifting
  # by whole turns puts it back at from, so p = 0 and p = 1 give from and from + 2 pi
  turns <- floor((from - q0(pf)) / (2 * pi) + 0.5)
  q + 2 * pi * turns
}

#' @rdname wrpcauchy
#' @export
#' @importFrom stats runif
rwrpcauchy <- function(n, mu = 0, rho, wrap = TRUE) {

  if (any(rho < 0) || any(rho >= 1)) stop("rho must be in the interval [0, 1).")

  n <- ceiling(n)
  angles <- qwrpcauchy(runif(n), mu = mu, rho = rho)

  # if generated angels should be wrapped, i.e. mapped to interval [-pi, pi], do so
  if(wrap){
    angles = (angles + pi) %% (2 * pi) - pi
  }
  angles
}
