#' Hurdle beta-binomial distribution
#'
#' Probability mass function and random generation for the hurdle (zero-altered) beta-binomial
#' distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' A hurdle distribution models the zeros and the positive counts as two separate
#' processes: the probability of a zero is a free parameter, and the positive counts
#' follow the corresponding zero-truncated distribution. Writing \eqn{p_0} for
#' \code{zeroprob},
#' \deqn{P(X = 0) = p_0, \qquad
#'       P(X = x) = (1 - p_0)\,\frac{P_{\mathrm{BB}}(x;\,n,a,b)}{1 - \pi_0},
#'       \quad x = 1, \ldots, n.}
#' where \eqn{\pi_0 = P_{\mathrm{BB}}(0;\,n,a,b)} is the probability of a zero under
#' the ordinary beta-binomial.
#'
#' Unlike zero-inflation, which can only add zeros, \code{zeroprob} here is exactly
#' the probability of observing a zero and may be larger \emph{or} smaller than the
#' beta-binomial would give on its own.
#'
#' Like \code{\link{betabinom}} itself, this distribution provides no distribution
#' function: the beta-binomial cdf has no closed form and would have to be summed
#' over the support, which cannot be taped for automatic differentiation.
#'
#' @seealso [betabinom], [zibetabinom], [ztbetabinom], [hbinom]
#'
#' @param x integer vector of counts
#' @param n number of random values to return.
#' @param size number of trials (zero or more)
#' @param shape1,shape2 positive shape parameters of the mixing beta distribution
#' @param zeroprob probability of a zero, between 0 and 1
#' @param log logical; return log-density if TRUE
#'
#' @return
#' \code{dhbetabinom} gives the probability mass function and \code{rhbetabinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rhbetabinom(5, size = 10, shape1 = 2, shape2 = 3, zeroprob = 0.4)
#' d <- dhbetabinom(x, size = 10, shape1 = 2, shape2 = 3, zeroprob = 0.4)
#' @name hbetabinom
NULL
#' @rdname hbetabinom
#' @export
dhbetabinom <- function(x, size, shape1, shape2, zeroprob = 0.5, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(size < 0) || any(size != floor(size))) stop("size must be non-negative integers.")
    if (any(shape1 <= 0) || any(shape2 <= 0)) stop("shape1 and shape2 must be positive.")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  if (inherits(x, "simref")) {
    return(dGenericSim("dhbetabinom", x = x, size = size, shape1 = shape1,
                       shape2 = shape2, zeroprob = zeroprob, log = log))
  }
  if (inherits(x, "osa")) {
    stop("Hurdle beta-binomial does not support OSA residuals.")
  }

  # log_hurdle() combines the point mass at zero with the rescaled zero-truncated
  # beta-binomial on the log scale, so either branch may be zero without a NaN
  # x is clamped because dbetabinom() hits the lgamma pole at negative integers;
  # log_hurdle()'s support indicators then remove the clamped values again
  logdens <- log_hurdle(x, dbetabinom(pmax.ad(x, 0), size, shape1, shape2, log = TRUE),
                        dbetabinom(0, size, shape1, shape2), zeroprob)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname hbetabinom
#' @export
#' @importFrom stats runif
rhbetabinom <- function(n, size, shape1, shape2, zeroprob = 0.5) {
  if (any(size < 0) || any(size != floor(size))) stop("size must be non-negative integers.")
  if (any(shape1 <= 0) || any(shape2 <= 0)) stop("shape1 and shape2 must be positive.")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  size <- rep_len(size, n); shape1 <- rep_len(shape1, n)
  shape2 <- rep_len(shape2, n); zeroprob <- rep_len(zeroprob, n)

  x <- numeric(n)
  cleared <- runif(n) >= zeroprob  # TRUE where the hurdle is cleared
  if (any(cleared)) x[cleared] <- rztbetabinom(sum(cleared), size[cleared],
                                               shape1[cleared], shape2[cleared])

  return(x)
}
