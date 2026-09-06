#' Johnson SU distribution (JSU)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Johnson SU distribution, in the original and in the moment parameterisation.
#'
#' @details
#' The Johnson SU distribution is a four-parameter continuous distribution on the
#' whole real line, obtained by applying the transformation
#' \deqn{Z = \nu + \tau \,\mathrm{asinh}\!\left(\frac{x-\mu}{\sigma}\right) \sim N(0,1)}
#' to a standard normal variable. It covers a wide range of skewness and kurtosis
#' combinations and is a common alternative to the Box-Cox families for data that
#' are not restricted to be positive.
#'
#' \code{djsu} uses the original parameterisation, in which \eqn{\mu} and
#' \eqn{\sigma} are a location and a scale parameter, \eqn{\nu} controls skewness
#' and \eqn{\tau} controls kurtosis. The density is
#' \deqn{f(x; \mu, \sigma, \nu, \tau) = \frac{\tau}{\sigma \sqrt{2\pi}}
#'   \frac{1}{\sqrt{z^2 + 1}} \exp\!\left(-\frac{r^2}{2}\right),}
#' where \eqn{z = (x-\mu)/\sigma} and \eqn{r = \nu + \tau\,\mathrm{asinh}(z)}.
#' Here \eqn{\mu} is \emph{not} the mean and \eqn{\sigma} is \emph{not} the standard deviation.
#'
#' \code{djsu2} uses the moment parameterisation, in which \eqn{\mu} \strong{is} the
#' mean and \eqn{\sigma} \strong{is} the standard deviation of the distribution, for
#' any admissible \eqn{\nu} and \eqn{\tau}. This is usually the more convenient
#' parameterisation for regression modelling, because the location and scale
#' parameters keep their interpretation as \eqn{\nu} and \eqn{\tau} change.
#' Writing \eqn{\omega = -\nu/\tau} and \eqn{w = \exp(\tau^{-2})}, it is obtained from
#' the original parameterisation by the reparameterisation
#' \deqn{c = \left[\tfrac{1}{2}(w-1)\left(w \cosh(2\omega) + 1\right)\right]^{-1/2},
#'   \qquad \sigma^* = c\,\sigma, \qquad
#'   \mu^* = \mu + c\,\sigma\sqrt{w}\,\sinh(\omega),}
#' so that \code{djsu2(x, mu, sigma, nu, tau)} equals
#' \code{djsu(x, }\eqn{\mu^*}\code{, }\eqn{\sigma^*}\code{, -nu, tau)}.
#'
#' These correspond to the \code{JSUo} and \code{JSU} families of the
#' \code{gamlss.dist} package respectively; see Chapter 18 of Rigby et al. (2019).
#' Note that the sign convention for \eqn{\nu} differs between the two
#' parameterisations, exactly as it does in \code{gamlss.dist}.
#'
#' All four \code{d} and \code{p} functions are compatible with automatic
#' differentiation by \code{RTMB}, so both simulation and one-step-ahead
#' residuals are supported.
#'
#' @references
#' Johnson, N. L. (1954). Systems of frequency curves derived from the first law of Laplace.
#' Trabajos de Estadistica, 5, 283-291.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [bcpe], [bct], [skewt], [skewnorm]
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter for \code{djsu}; the mean for \code{djsu2}.
#' @param sigma scale parameter for \code{djsu}; the standard deviation for \code{djsu2}. Must be positive.
#' @param nu skewness parameter (real). Positive \eqn{\nu} gives left skewness in \code{djsu} and right skewness in \code{djsu2}.
#' @param tau kurtosis parameter, must be positive. Large \eqn{\tau} approaches the normal distribution.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{djsu} gives the density, \code{pjsu} gives the distribution function, \code{qjsu} gives the quantile function, and \code{rjsu} generates random deviates.
#' \code{djsu2}, \code{pjsu2}, \code{qjsu2} and \code{rjsu2} are the corresponding functions for the moment parameterisation.
#'
#' @examples
#' set.seed(123)
#' # original parameterisation
#' x <- rjsu(5, mu = 0, sigma = 1, nu = -1, tau = 2)
#' d <- djsu(x, mu = 0, sigma = 1, nu = -1, tau = 2)
#' p <- pjsu(x, mu = 0, sigma = 1, nu = -1, tau = 2)
#' q <- qjsu(p, mu = 0, sigma = 1, nu = -1, tau = 2)
#'
#' # moment parameterisation: mu is the mean, sigma the standard deviation
#' y <- rjsu2(1000, mu = 3, sigma = 2, nu = 1, tau = 3)
#' c(mean = mean(y), sd = sd(y))
#' @name jsu
NULL

#' @rdname jsu
#' @export
#' @import RTMB
djsu <- function(x, mu = 0, sigma = 1, nu = 0, tau = 1, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(tau <= 0)) stop("tau must be > 0")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("djsu", x = x, mu = mu, sigma = sigma, nu = nu, tau = tau, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("djsu", x = x, mu = mu, sigma = sigma, nu = nu, tau = tau, log = log))
  }

  z <- (x - mu) / sigma
  r <- nu + tau * asinh(z)

  logdens <- log(tau) - log(sigma) - 0.5 * log1p(z * z) - 0.5 * log(2 * pi) - 0.5 * r * r

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname jsu
#' @export
pjsu <- function(q, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(tau <= 0)) stop("tau must be > 0")
  }

  # the transformation to normality is exact, so the cdf is the normal cdf of r
  r <- nu + tau * asinh((q - mu) / sigma)
  p <- RTMB::pnorm(r)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}

#' @rdname jsu
#' @export
qjsu <- function(p, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(tau <= 0)) stop("tau must be > 0")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  mu + sigma * sinh((stats::qnorm(p) - nu) / tau)
}

#' @rdname jsu
#' @export
#' @importFrom stats runif
rjsu <- function(n, mu = 0, sigma = 1, nu = 0, tau = 1) {

  if (any(sigma <= 0)) stop("sigma must be > 0")
  if (any(tau <= 0)) stop("tau must be > 0")

  n <- ceiling(n)
  p <- runif(n)

  qjsu(p, mu = mu, sigma = sigma, nu = nu, tau = tau)
}

# internal helper: moment parameterisation (mean, sd) -> original (location, scale).
# expm1() is used for w - 1 because for large tau the argument tau^-2 is tiny and
# exp(tau^-2) - 1 cancels; gamlss.dist instead switches branch at tau > 1e7.
.jsu2_internal_params <- function(mu, sigma, nu, tau) {

  rtau <- 1 / tau
  w <- exp(rtau * rtau)
  wm1 <- expm1(rtau * rtau)
  omega <- -nu * rtau

  cc <- 1 / sqrt(0.5 * wm1 * (w * cosh(2 * omega) + 1))

  list(location = mu + cc * sigma * sqrt(w) * sinh(omega),
       scale    = cc * sigma)
}

#' @rdname jsu
#' @export
djsu2 <- function(x, mu = 0, sigma = 1, nu = 0, tau = 1, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(tau <= 0)) stop("tau must be > 0")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("djsu2", x = x, mu = mu, sigma = sigma, nu = nu, tau = tau, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("djsu2", x = x, mu = mu, sigma = sigma, nu = nu, tau = tau, log = log))
  }

  pars <- .jsu2_internal_params(mu, sigma, nu, tau)

  djsu(x, mu = pars$location, sigma = pars$scale, nu = -nu, tau = tau, log = log)
}

#' @rdname jsu
#' @export
pjsu2 <- function(q, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(tau <= 0)) stop("tau must be > 0")
  }

  pars <- .jsu2_internal_params(mu, sigma, nu, tau)

  pjsu(q, mu = pars$location, sigma = pars$scale, nu = -nu, tau = tau,
       lower.tail = lower.tail, log.p = log.p)
}

#' @rdname jsu
#' @export
qjsu2 <- function(p, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(tau <= 0)) stop("tau must be > 0")
  }

  pars <- .jsu2_internal_params(mu, sigma, nu, tau)

  qjsu(p, mu = pars$location, sigma = pars$scale, nu = -nu, tau = tau,
       lower.tail = lower.tail, log.p = log.p)
}

#' @rdname jsu
#' @export
#' @importFrom stats runif
rjsu2 <- function(n, mu = 0, sigma = 1, nu = 0, tau = 1) {

  if (any(sigma <= 0)) stop("sigma must be > 0")
  if (any(tau <= 0)) stop("tau must be > 0")

  n <- ceiling(n)
  p <- runif(n)

  qjsu2(p, mu = mu, sigma = sigma, nu = nu, tau = tau)
}
