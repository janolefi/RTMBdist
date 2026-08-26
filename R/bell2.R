#' Reparameterised Bell distribution
#'
#' Probability mass function, distribution function, quantile function, and
#' random generation for the Bell distribution reparameterised in terms of its
#' mean.
#'
#' @details
#' This implementation of \code{dbell2} and \code{pbell2} allows for automatic
#' differentiation with \code{RTMB} with respect to \code{mu}.
#'
#' The Bell distribution has mean \eqn{\mu = \theta e^{\theta}}, which is a
#' bijection from \eqn{\theta > 0} to \eqn{\mu > 0} and is inverted by the
#' principal branch of the Lambert W function,
#' \deqn{\theta = W(\mu).}
#' Every positive mean therefore corresponds to exactly one \eqn{\theta}. All
#' four functions simply apply this transformation and hand over to their
#' \code{\link{bell}} counterparts.
#'
#' In this parameterisation the variance is \eqn{\mu (1 + W(\mu))}, so the
#' distribution is always overdispersed relative to the Poisson distribution,
#' but the degree of overdispersion is determined by the mean rather than by a
#' free parameter.
#'
#' \code{\link{lambertW}} is AD-compatible to arbitrary order, so \code{mu} may
#' be a parameter of a model fitted by Laplace approximation.
#'
#' @param x,q integer vector of counts
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mu vector of positive means
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dbell2} gives the probability mass function, \code{pbell2} gives the distribution function, \code{qbell2} gives the quantile function, and \code{rbell2} generates random deviates.
#'
#' @references
#' Castellares, F., Ferrari, S. L. P., and Lemonte, A. J. (2018). On the Bell
#' distribution and its associated regression model for count data.
#' \emph{Applied Mathematical Modelling} 56, 172-185.
#' doi:10.1016/j.apm.2017.12.014
#'
#' @seealso \code{\link{bell}} for the natural parameterisation.
#'
#' @examples
#' set.seed(123)
#' x <- rbell2(1, 3)
#' d <- dbell2(x, 3)
#' p <- pbell2(x, 3)
#' q <- qbell2(p, 3)
#'
#' # the two parameterisations agree
#' all.equal(dbell2(0:5, 3), dbell(0:5, lambertW(3)))
#' @name bell2
NULL

#' @rdname bell2
#' @export
#' @import RTMB
dbell2 <- function(x, mu, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(mu <= 0)) stop("mu must be > 0")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dbell2", x = x, mu = mu, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dbell2", x = x, mu = mu, log = log))
  }

  dbell(x, theta = lambertW(mu), log = log)
}

#' @rdname bell2
#' @export
#' @import RTMB
pbell2 <- function(q, mu, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(mu <= 0)) stop("mu must be > 0")
  }

  pbell(q, theta = lambertW(mu), lower.tail = lower.tail, log.p = log.p)
}

#' @rdname bell2
#' @export
qbell2 <- function(p, mu, lower.tail = TRUE, log.p = FALSE) {

  if (any(mu <= 0)) stop("mu must be > 0")

  qbell(p, theta = lambertW(mu), lower.tail = lower.tail, log.p = log.p)
}

#' @rdname bell2
#' @export
rbell2 <- function(n, mu) {

  if (any(mu <= 0)) stop("mu must be > 0")

  rbell(n, theta = lambertW(mu))
}
