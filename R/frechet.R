#' Frechet distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the Frechet distribution.
#'
#' @details
#' \code{dfrechet} and \code{pfrechet} allow for automatic differentiation with \code{RTMB}.
#'
#' With \eqn{y = (x - \mu) / \sigma} the density is
#' \deqn{f(x;\,\mu,\sigma,\alpha) = \frac{\alpha}{\sigma} y^{-\alpha - 1} \exp(-y^{-\alpha}), \quad x > \mu,}
#' and the distribution function is \eqn{F(x) = \exp(-y^{-\alpha})}.
#'
#' The Frechet distribution is the heavy-tailed extreme value distribution: it
#' is the \link[=gev]{generalised extreme value} distribution with shape
#' \eqn{\xi = 1/\alpha > 0}, reparameterised so that the shape enters as a tail
#' index rather than as a reciprocal. All moments of order \eqn{\alpha} and
#' above are infinite, so the mean exists only for \eqn{\alpha > 1} and the
#' variance only for \eqn{\alpha > 2}.
#'
#' @references
#' Frechet, M. (1927) Sur la loi de probabilite de l'ecart maximum. Annales de
#' la Societe Polonaise de Mathematique, 6, 93-116.
#'
#' Kotz, S. and Nadarajah, S. (2000) Extreme Value Distributions: Theory and
#' Applications, Imperial College Press, doi:10.1142/p191.
#'
#' @seealso [gev], [gumbel], [pareto]
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter, the lower end point of the support.
#' @param sigma scale parameter, must be positive.
#' @param alpha shape parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dfrechet} gives the density, \code{pfrechet} gives the distribution function, \code{qfrechet} gives the quantile function, and \code{rfrechet} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rfrechet(5, mu = 0, sigma = 1, alpha = 3)
#' d <- dfrechet(x, mu = 0, sigma = 1, alpha = 3)
#' p <- pfrechet(x, mu = 0, sigma = 1, alpha = 3)
#' q <- qfrechet(p, mu = 0, sigma = 1, alpha = 3)
#' @name frechet
NULL

#' @rdname frechet
#' @export
dfrechet <- function(x, mu = 0, sigma = 1, alpha = 1, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if(any(sigma <= 0)) stop("sigma must be > 0")
    if(any(alpha <= 0)) stop("alpha must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dfrechet", x=x, mu=mu, sigma=sigma, alpha=alpha, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dfrechet", x=x, mu=mu, sigma=sigma, alpha=alpha, log=log))
  }

  y <- (x - mu) / sigma
  insup <- greater(y, 0) # the support is x > mu

  # outside the support y is replaced by 1, so that log y and y^-alpha stay
  # finite; the density itself is set to zero by log(insup)
  logy <- log(insup * y + (1 - insup))

  logdens <- log(insup) + log(alpha) - log(sigma) -
    (alpha + 1) * logy - exp(-alpha * logy)

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname frechet
#' @export
pfrechet <- function(q, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if(any(sigma <= 0)) stop("sigma must be > 0")
    if(any(alpha <= 0)) stop("alpha must be > 0")
  }

  y <- (q - mu) / sigma
  insup <- greater(y, 0)

  logy <- log(insup * y + (1 - insup))

  p <- insup * exp(-exp(-alpha * logy))

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  return(p)
}

#' @rdname frechet
#' @export
qfrechet <- function(p, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if(any(sigma <= 0)) stop("sigma must be > 0")
    if(any(alpha <= 0)) stop("alpha must be > 0")
  }

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p

  if(!ad_context()) {
    if(any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  mu + sigma * (-log(p))^(-1 / alpha)
}

#' @rdname frechet
#' @export
#' @importFrom stats runif
rfrechet <- function(n, mu = 0, sigma = 1, alpha = 1) {

  if(!ad_context()) {
    if(any(sigma <= 0)) stop("sigma must be > 0")
    if(any(alpha <= 0)) stop("alpha must be > 0")
  }

  n <- ceiling(n)
  p <- runif(n)

  qfrechet(p, mu = mu, sigma = sigma, alpha = alpha)
}
