#' Hurdle negative binomial distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the hurdle (zero-altered) negative binomial distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' A hurdle distribution models the zeros and the positive counts as two separate
#' processes: the probability of a zero is a free parameter, and the positive
#' counts follow the corresponding zero-truncated distribution. Writing
#' \eqn{p_0} for \code{zeroprob},
#' \deqn{P(X = 0) = p_0, \qquad
#'       P(X = x) = (1 - p_0)\,\frac{P_{\mathrm{NB}}(x;\,r,\pi)}{1 - \pi^{r}},
#'       \quad x = 1, 2, \ldots}
#'
#' Unlike zero-inflation, which can only add zeros to those the negative binomial
#' already produces, \code{zeroprob} here is exactly the probability of a zero and
#' may be larger \emph{or} smaller than \eqn{\pi^{r}}. The two coincide with the
#' ordinary negative binomial when \code{zeroprob} equals \eqn{\pi^{r}}.
#'
#' @references
#' Mullahy, J. (1986) Specification and testing of some modified count data models.
#' Journal of Econometrics, 33, 341-365.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [hnbinom2], [hpois], [hbinom], [zinbinom], [ztnbinom]
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param size dispersion parameter, must be strictly positive
#' @param prob probability of success in each trial, in (0,1)
#' @param zeroprob probability of a zero, between 0 and 1
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dhnbinom} gives the probability mass function, \code{phnbinom} gives the distribution function, and \code{rhnbinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rhnbinom(5, size = 2, prob = 0.4, zeroprob = 0.3)
#' d <- dhnbinom(x, size = 2, prob = 0.4, zeroprob = 0.3)
#' p <- phnbinom(x, size = 2, prob = 0.4, zeroprob = 0.3)
#' @name hnbinom
NULL
#' @rdname hnbinom
#' @export
#' @importFrom RTMB dnbinom logspace_add
dhnbinom <- function(x, size, prob, zeroprob = 0.5, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(size <= 0)) stop("size must be > 0")
    if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dhnbinom", x = x, size = size, prob = prob, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dhnbinom", x = x, size = size, prob = prob, zeroprob = zeroprob, log=log))
  }

  # the point mass at zero and the rescaled zero-truncated negative binomial are
  # combined on the log scale, so either branch may be exactly zero without a NaN
  log_zero <- log(zeroprob) + log(iszero(x))
  log_pos <- log1p(-zeroprob) + dnbinom(x, size = size, prob = prob, log = TRUE) -
    log1p(-dnbinom(0, size = size, prob = prob)) + log(ispos_strict(x))

  logdens <- logspace_add(log_zero, log_pos)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname hnbinom
#' @export
phnbinom <- function(q, size, prob, zeroprob = 0.5, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if (any(size <= 0)) stop("size must be > 0")
    if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    q <- floor(q)  # make sure it's integer-valued
  }

  # zero mass plus the rescaled zero-truncated cdf, and zero below the support
  p <- greater(q, -1) * (zeroprob + (1 - zeroprob) * pztnbinom(q, size, prob))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname hnbinom
#' @export
#' @importFrom stats runif
rhnbinom <- function(n, size, prob, zeroprob = 0.5) {
  if (any(size <= 0)) stop("size must be > 0")
  if (any(prob <= 0 | prob >= 1)) stop("prob must be in (0,1)")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  size <- rep_len(size, n)
  prob <- rep_len(prob, n)
  zeroprob <- rep_len(zeroprob, n)

  x <- numeric(n)
  cleared <- runif(n) >= zeroprob  # TRUE where the hurdle is cleared
  if (any(cleared)) x[cleared] <- rztnbinom(sum(cleared), size[cleared], prob[cleared])

  return(x)
}
