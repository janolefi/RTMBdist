#' Zero-inflated and reparameterised gamma distribution
#'
#' Density, distribution function, and random generation for
#' the zero-inflated gamma distribution reparameterised in terms of mean and standard deviation.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' Uses the same density as \code{zigamma} with \eqn{\text{shape} = \mu^2/s^2} and \eqn{\text{scale} = s^2/\mu}:
#' \deqn{f(x;\,\mu,s,p_0) = p_0\,\mathbf{1}[x=0] + (1-p_0)\,f_{\mathrm{Gamma}}(x;\,\mu^2/s^2,\,s^2/\mu)\,\mathbf{1}[x>0].}
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return
#' @param mean mean parameter, must be positive.
#' @param sd standard deviation parameter, must be positive.
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{dzigamma2} gives the density, \code{pzigamma2} gives the distribution function, and \code{rzigamma2} generates random deviates.
#'
#' @examples
#' x <- rzigamma2(1, 2, 1, 0.5)
#' d <- dzigamma2(x, 2, 1, 0.5)
#' p <- pzigamma2(x, 2, 1, 0.5)
#' @name zigamma2
NULL

#' @rdname zigamma2
#' @export
#' @importFrom RTMB dgamma logspace_add
dzigamma2 = function(x, mean = 1, sd = 1, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure shape > 0, scale > 0, zeroprob in [0,1]
    if (any(mean <= 0)) stop("mean must be > 0")
    if (any(sd <= 0)) stop("sd must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # parameter transformation
  shape = mean^2 / sd^2
  scale = sd^2 / mean

  dzigamma(x, shape = shape, scale = scale, zeroprob = zeroprob, log = log)
}
#' @rdname zigamma2
#' @export
pzigamma2 <- function(q, mean = 1, sd = 1, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {
  if(!ad_context()) {
    # ensure shape > 0, scale > 0, zeroprob in [0,1]
    if (any(mean <= 0)) stop("mean must be > 0")
    if (any(sd <= 0)) stop("sd must be > 0")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # parameter transformation
  shape = mean^2 / sd^2
  scale = sd^2 / mean

  pzigamma(q, shape = shape, scale = scale, zeroprob = zeroprob,
           lower.tail = lower.tail, log.p = log.p)
}
#' @rdname zigamma2
#' @export
rzigamma2 <- function(n, mean = 1, sd = 1, zeroprob = 0) {
  # ensure mean, sd > 0
  if (any(mean <= 0)) stop("mean must be strictly positive.")
  if (any(sd <= 0)) stop("sd must be strictly positive.")
  # ensure zeroprob in [0,1]
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  # parameter transformation
  shape = mean^2 / sd^2
  scale = sd^2 / mean

  rzigamma(n, shape = shape, scale = scale, zeroprob = zeroprob)
}
