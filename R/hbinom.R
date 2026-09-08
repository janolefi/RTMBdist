#' Hurdle binomial distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the hurdle (zero-altered) binomial distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' A hurdle distribution models the zeros and the positive counts as two separate
#' processes: the probability of a zero is a free parameter, and the positive
#' counts follow the corresponding zero-truncated distribution. Writing
#' \eqn{p_0} for \code{zeroprob},
#' \deqn{P(X = 0) = p_0, \qquad
#'       P(X = x) = (1 - p_0)\,\frac{P_{\mathrm{Bin}}(x;\,n,\pi)}{1 - \pi_0},
#'       \quad x = 1, \ldots, n.}
#' where \eqn{\pi_0 = P_{\mathrm{Bin}}(0;\,n,\pi)} is the probability of a zero under
#' the ordinary binomial.
#'
#' Unlike zero-inflation, which can only add zeros to those the binomial already
#' produces, \code{zeroprob} here is exactly the probability of a zero and may be
#' larger \emph{or} smaller than \eqn{\pi_0}. The two coincide with the
#' ordinary binomial when \code{zeroprob} equals \eqn{\pi_0}.
#'
#' @references
#' Mullahy, J. (1986) Specification and testing of some modified count data models.
#' Journal of Econometrics, 33, 341-365.
#'
#' Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC,
#' doi:10.1201/9780429298547. An older version can be found in https://www.gamlss.com/.
#'
#' @seealso [hpois], [hnbinom], [zibinom], [ztbinom]
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param size number of trials (zero or more)
#' @param prob probability of success on each trial
#' @param zeroprob probability of a zero, between 0 and 1
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dhbinom} gives the probability mass function, \code{phbinom} gives the distribution function, and \code{rhbinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rhbinom(5, size = 10, prob = 0.3, zeroprob = 0.4)
#' d <- dhbinom(x, size = 10, prob = 0.3, zeroprob = 0.4)
#' p <- phbinom(x, size = 10, prob = 0.3, zeroprob = 0.4)
#' @name hbinom
NULL
#' @rdname hbinom
#' @export
#' @importFrom RTMB dbinom logspace_add
dhbinom <- function(x, size, prob, zeroprob = 0.5, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    if (any(size < 0)) stop("size must be >= 0")
    if (any(prob < 0 | prob > 1)) stop("prob must be in [0,1]")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dhbinom", x = x, size = size, prob = prob, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dhbinom", x = x, size = size, prob = prob, zeroprob = zeroprob, log=log))
  }

  # log_hurdle() combines the point mass at zero with the rescaled zero-truncated
  # binomial on the log scale, so either branch may be exactly zero without a NaN
  logdens <- log_hurdle(x, dbinom(x, size, prob, log = TRUE),
                        dbinom(0, size, prob), zeroprob)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname hbinom
#' @export
phbinom <- function(q, size, prob, zeroprob = 0.5, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if (any(size < 0)) stop("size must be >= 0")
    if (any(prob < 0 | prob > 1)) stop("prob must be in [0,1]")
    if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
    q <- floor(q)  # make sure it's integer-valued
  }

  # zero mass plus the rescaled zero-truncated cdf, and zero below the support
  p <- greater(q, -1) * (zeroprob + (1 - zeroprob) * pztbinom(q, size, prob))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname hbinom
#' @export
#' @importFrom stats runif
rhbinom <- function(n, size, prob, zeroprob = 0.5) {
  if (any(size < 0)) stop("size must be >= 0")
  if (any(prob < 0 | prob > 1)) stop("prob must be in [0,1]")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  size <- rep_len(size, n)
  prob <- rep_len(prob, n)
  zeroprob <- rep_len(zeroprob, n)

  x <- numeric(n)
  cleared <- runif(n) >= zeroprob  # TRUE where the hurdle is cleared
  if (any(cleared)) x[cleared] <- rztbinom(sum(cleared), size[cleared], prob[cleared])

  return(x)
}
