#' Reparameterised hurdle negative binomial distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the hurdle (zero-altered) negative binomial distribution, parameterised by the
#' mean of the untruncated negative binomial.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' This is \code{\link{hnbinom}} with the success probability replaced by
#' \deqn{\pi = \frac{\mathrm{size}}{\mathrm{size} + \mu},}
#' so that \eqn{\mu} is the mean of the \emph{untruncated} negative binomial, whose
#' variance is \eqn{\mu + \mu^2/\mathrm{size}}. Note that \eqn{\mu} is not the mean
#' of the hurdle distribution itself, which also depends on \code{zeroprob}.
#'
#' As for all hurdle distributions, \code{zeroprob} is exactly the probability of
#' observing a zero and may be larger \emph{or} smaller than the negative binomial
#' would give on its own.
#'
#' @references
#' Mullahy, J. (1986) Specification and testing of some modified count data models.
#' Journal of Econometrics, 33, 341-365.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [hnbinom], [hpois], [nbinom2], [zinbinom2], [ztnbinom2]
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param mu mean of the untruncated negative binomial, must be strictly positive
#' @param size dispersion parameter, must be strictly positive
#' @param zeroprob probability of a zero, between 0 and 1
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dhnbinom2} gives the probability mass function, \code{phnbinom2} gives the distribution function, and \code{rhnbinom2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rhnbinom2(5, mu = 3, size = 2, zeroprob = 0.3)
#' d <- dhnbinom2(x, mu = 3, size = 2, zeroprob = 0.3)
#' p <- phnbinom2(x, mu = 3, size = 2, zeroprob = 0.3)
#' @name hnbinom2
NULL
#' @rdname hnbinom2
#' @export
dhnbinom2 <- function(x, mu, size, zeroprob = 0.5, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(size <= 0)) stop("size must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dhnbinom2", x = x, mu = mu, size = size, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dhnbinom2", x = x, mu = mu, size = size, zeroprob = zeroprob, log=log))
  }

  prob <- size / (size + mu)

  dhnbinom(x, size = size, prob = prob, zeroprob = zeroprob, log = log)
}
#' @rdname hnbinom2
#' @export
phnbinom2 <- function(q, mu, size, zeroprob = 0.5, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(size <= 0)) stop("size must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    q <- floor(q)  # make sure it's integer-valued
  }

  prob <- size / (size + mu)

  phnbinom(q, size = size, prob = prob, zeroprob = zeroprob,
           lower.tail = lower.tail, log.p = log.p)
}
#' @rdname hnbinom2
#' @export
rhnbinom2 <- function(n, mu, size, zeroprob = 0.5) {
  if (any(mu <= 0)) stop("mu must be > 0")
  if (any(size <= 0)) stop("size must be > 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  prob <- size / (size + mu)

  rhnbinom(n, size = size, prob = prob, zeroprob = zeroprob)
}
