#' Hurdle geometric distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the hurdle (zero-altered) geometric distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' A hurdle distribution models the zeros and the positive counts as two separate
#' processes: the probability of a zero is a free parameter, and the positive
#' counts follow the corresponding zero-truncated distribution. Writing
#' \eqn{p_0} for \code{zeroprob},
#' \deqn{P(X = 0) = p_0, \qquad
#'       P(X = x) = (1 - p_0)\,\frac{P_{\mathrm{Geom}}(x;\,\pi)}{1 - \pi_0},
#'       \quad x = 1, 2, \ldots}
#' where \eqn{\pi_0 = P_{\mathrm{Geom}}(0;\,\pi) = \pi} is the probability of a zero
#' under the ordinary geometric.
#'
#' Unlike zero-inflation, which can only add zeros to those the geometric already
#' produces, \code{zeroprob} here is exactly the probability of a zero and may be
#' larger \emph{or} smaller than \eqn{\pi_0}. The two coincide with the ordinary
#' geometric when \code{zeroprob} equals \eqn{\pi_0}.
#'
#' @references
#' Mullahy, J. (1986) Specification and testing of some modified count data models.
#' Journal of Econometrics, 33, 341-365.
#'
#' @seealso [zigeom], [ztgeom], [hnbinom], [hpois]
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param prob probability of success in each trial, in (0,1)
#' @param zeroprob probability of a zero, between 0 and 1
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dhgeom} gives the probability mass function, \code{phgeom} gives the distribution function, and \code{rhgeom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rhgeom(5, prob = 0.3, zeroprob = 0.4)
#' d <- dhgeom(x, prob = 0.3, zeroprob = 0.4)
#' p <- phgeom(x, prob = 0.3, zeroprob = 0.4)
#' @name hgeom
NULL
#' @rdname hgeom
#' @export
#' @import RTMB
dhgeom <- function(x, prob, zeroprob = 0.5, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dhgeom", x = x, prob = prob, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dhgeom", x = x, prob = prob, zeroprob = zeroprob, log=log))
  }

  # log_hurdle() combines the point mass at zero with the rescaled zero-truncated
  # geometric on the log scale, so either branch may be exactly zero without a NaN
  # log_hurdle() combines the point mass at zero with the rescaled zero-truncated
  # geometric on the log scale; P(X = 0) for the geometric is prob
  logdens <- log_hurdle(x, dgeom(x, prob, log = TRUE), prob, zeroprob)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname hgeom
#' @export
phgeom <- function(q, prob, zeroprob = 0.5, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    q <- floor(q)  # make sure it's integer-valued
  }

  # zero mass plus the rescaled zero-truncated cdf, and zero below the support
  p <- greater(q, -1) * (zeroprob + (1 - zeroprob) * pztgeom(q, prob))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname hgeom
#' @export
#' @importFrom stats runif
rhgeom <- function(n, prob, zeroprob = 0.5) {
  if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  prob <- rep_len(prob, n)
  zeroprob <- rep_len(zeroprob, n)

  x <- numeric(n)
  cleared <- runif(n) >= zeroprob  # TRUE where the hurdle is cleared
  if (any(cleared)) x[cleared] <- rztgeom(sum(cleared), prob[cleared])

  return(x)
}
