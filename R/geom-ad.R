#' AD-compatible geometric distribution
#'
#' Density and distribution function for the geometric distribution, written so
#' that they can be taped by \code{RTMB}.
#'
#' @details
#' \code{stats} already provides the geometric distribution, but its versions
#' cannot be differentiated. These are AD-compatible replacements, reached
#' automatically whenever an argument is an AD variable, so \code{stats::dgeom}
#' and \code{stats::pgeom} are left untouched for ordinary use.
#'
#' The parameterisation is the same as in \code{stats}: \eqn{X} is the number of
#' failures before the first success, so
#' \deqn{P(X = x;\,\pi) = \pi\,(1 - \pi)^{x}, \quad x = 0, 1, 2, \ldots}
#' The density is obtained from the negative binomial with \code{size = 1}, for
#' which \code{RTMB} provides an AD method; the distribution function
#' \eqn{1 - (1-\pi)^{x+1}} is elementary.
#'
#' @seealso [zigeom], [ztgeom], [hgeom]
#'
#' @param x,q integer vector of counts
#' @param prob probability of success in each trial, in (0,1]
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dgeom.ad} gives the probability mass function and \code{pgeom.ad} gives the distribution function.
#'
#' @examples
#' dgeom.ad(0:5, prob = 0.3)
#' pgeom.ad(0:5, prob = 0.3)
#' @name geom.ad
NULL

#' @rdname geom.ad
#' @export
#' @import RTMB
dgeom.ad <- function(x, prob, log = FALSE) {

  # the geometric is the negative binomial with size = 1, for which RTMB
  # provides an AD-compatible S4 method
  dnbinom(x, size = 1, prob = prob, log = log)
}

#' @rdname geom.ad
#' @export
pgeom.ad <- function(q, prob, lower.tail = TRUE, log.p = FALSE) {

  # P(X <= q) = 1 - (1 - prob)^(q + 1), and zero below the support
  p <- greater(q, -1) * (1 - (1 - prob)^(q + 1))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
