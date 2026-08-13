#' Conway-Maxwell-binomial distribution
#'
#' Probability mass function, distribution function, quantile function, and
#' random generation for the Conway-Maxwell-binomial (CMB) distribution.
#'
#' @details
#' This implementation of \code{dcombinom} and \code{pcombinom} allows for
#' automatic differentiation with \code{RTMB}, including differentiation with
#' respect to \code{x}, such that one-step-ahead (OSA) residuals are supported.
#'
#' The CMB distribution generalises the binomial distribution by an additional
#' dispersion parameter \eqn{\nu} in the same way that the
#' Conway-Maxwell-Poisson distribution generalises the Poisson distribution.
#' Its probability mass function is
#'
#' \deqn{P(X = x;\, n, p, \nu) = \frac{1}{Z(n, p, \nu)} \binom{n}{x}^{\nu} p^x (1-p)^{n-x}, \quad x = 0, 1, \ldots, n,}
#'
#' with normalising constant
#'
#' \deqn{Z(n, p, \nu) = \sum_{k=0}^{n} \binom{n}{k}^{\nu} p^k (1-p)^{n-k}.}
#'
#' For \eqn{\nu = 1} this reduces to the binomial distribution. Values
#' \eqn{\nu > 1} give under-dispersion and \eqn{\nu < 1} over-dispersion
#' relative to a binomial distribution with the same mean. As the support is
#' finite, \eqn{Z} converges for every real \eqn{\nu}, so \eqn{\nu} is not
#' restricted to be positive; negative values yield strongly over-dispersed,
#' U-shaped distributions.
#'
#' The distribution arises as the sum of \eqn{n} exchangeable, possibly
#' associated Bernoulli variables, where \eqn{\nu} controls the association.
#' Note that \code{prob} is \emph{not} the mean divided by \code{size} unless
#' \eqn{\nu = 1}; the mean has no closed form for general \eqn{\nu} and must be
#' obtained by summation over the support.
#'
#' @param x,q integer vector of counts in \eqn{\{0, 1, \ldots, }\code{size}\eqn{\}}
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param size vector of numbers of trials (non-negative integers)
#' @param prob vector of success probabilities in \eqn{(0, 1)}
#' @param nu vector of dispersion parameters; \code{nu = 1} gives the binomial
#'   distribution, \code{nu > 1} under-dispersion and \code{nu < 1}
#'   over-dispersion. May be negative.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities are
#'   returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dcombinom} gives the probability mass function, \code{pcombinom} gives
#' the distribution function, \code{qcombinom} gives the quantile function, and
#' \code{rcombinom} generates random deviates.
#'
#' @references
#' Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S., and Boatwright, P.
#' (2005). A useful distribution for fitting discrete data: revival of the
#' Conway-Maxwell-Poisson distribution. \emph{Journal of the Royal Statistical
#' Society: Series C} 54(1), 127-142.
#'
#' Kadane, J. B. (2016). Sums of possibly associated Bernoulli variables: the
#' Conway-Maxwell-binomial distribution. \emph{Bayesian Analysis} 11(2),
#' 403-420.
#'
#' \url{https://en.wikipedia.org/wiki/Conway-Maxwell-binomial_distribution}
#'
#' @examples
#' set.seed(123)
#' x <- rcombinom(1, size = 10, prob = 0.4, nu = 1.5)
#' d <- dcombinom(x, 10, 0.4, 1.5)
#' p <- pcombinom(x, 10, 0.4, 1.5)
#' q <- qcombinom(p, 10, 0.4, 1.5)
#'
#' # nu = 1 recovers the binomial distribution
#' all.equal(dcombinom(0:10, 10, 0.3, 1), dbinom(0:10, 10, 0.3))
#' @name combinom
NULL

# Recycle to length n, using integer indexing so that advectors are handled.
recycle_combinom <- function(x, n) {
  if (length(x) == n) return(x)
  x[rep_len(seq_along(x), n)]
}

# log C(size, x) via lgamma, so that it is smooth (and AD-able) in x.
lchoose_combinom <- function(size, x) {
  lgamma(size + 1) - lgamma(x + 1) - lgamma(size - x + 1)
}

# log normalising constant log Z, folded with logspace_add for stability.
# size must be a single numeric value; psi and nu are vectors of equal length.
lZ_combinom <- function(psi, nu, size) {
  lchoosek <- lchoose(size, 0:size) # numeric constants: size is always data
  out <- nu * lchoosek[1] + 0 * psi # k = 0 term, log C(size, 0) = 0
  for (k in seq_len(size)) {
    out <- RTMB::logspace_add(out, nu * lchoosek[k + 1] + k * psi)
  }
  out
}

# Shared argument checking for all four functions.
check_combinom_pars <- function(size, prob, nu) {
  if (any(size < 0) || any(size != floor(size)))
    stop("size must be non-negative integers.")
  if (any(prob <= 0) || any(prob >= 1))
    stop("prob must be in (0, 1).")
  if (any(!is.finite(nu)))
    stop("nu must be finite.")
  invisible(NULL)
}

#' @rdname combinom
#' @export
#' @import RTMB
dcombinom <- function(x, size, prob, nu = 1, log = FALSE) {

  if (inherits(size, "advector"))
    stop("size must be numeric data, not a parameter.")

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    check_combinom_pars(size, prob, nu)
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dcombinom", x = x, size = size, prob = prob, nu = nu, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dcombinom", x = x, size = size, prob = prob, nu = nu, log = log))
  }

  nx <- max(length(x), length(size), length(prob), length(nu))
  x <- recycle_combinom(x, nx)
  size <- recycle_combinom(size, nx)
  prob <- recycle_combinom(prob, nx)
  nu <- recycle_combinom(nu, nx)

  psi <- log(prob) - log1p(-prob)

  # log Z depends on size, so evaluate it once per distinct size. The
  # accumulator must involve both psi and nu so that it is an advector
  # whenever either of them is one.
  lZ <- 0 * psi + 0 * nu
  for (s in unique(size)) {
    idx <- which(size == s)
    lZ[idx] <- lZ_combinom(psi[idx], nu[idx], s)
  }

  logdens <- nu * lchoose_combinom(size, x) + x * psi - lZ

  # The lgamma form above extends the pmf smoothly to non-integer x, which is
  # what makes it differentiable in x for OSA residuals. Outside of AD the
  # density of a non-integer is zero; x < 0 and x > size already give -Inf.
  if (!inherits(x, "advector")) {
    noninteger <- (x != floor(x))
    if (any(noninteger)) logdens[noninteger] <- -Inf
  }

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname combinom
#' @export
#' @import RTMB
pcombinom <- function(q, size, prob, nu = 1, lower.tail = TRUE, log.p = FALSE) {

  if (inherits(size, "advector"))
    stop("size must be numeric data, not a parameter.")

  if (!ad_context()) {
    check_combinom_pars(size, prob, nu)
  }

  nq <- max(length(q), length(size), length(prob), length(nu))
  q <- recycle_combinom(q, nq)
  size <- recycle_combinom(size, nq)
  prob <- recycle_combinom(prob, nq)
  nu <- recycle_combinom(nu, nq)

  psi <- log(prob) - log1p(-prob)

  # The CDF is a sum of pmf values over the support, masked by the smooth
  # indicator 1{k <= q}. Writing it this way (rather than indexing a cumulative
  # sum) keeps it evaluable when q is an advector, as required for OSA.
  # The half-integer shift avoids the ambiguity of ispos() at exactly zero.
  p <- 0 * q + 0 * psi + 0 * nu
  for (s in unique(size)) {
    idx <- which(size == s)
    lZ <- lZ_combinom(psi[idx], nu[idx], s)
    ps <- 0 * lZ
    for (k in 0:s) {
      ps <- ps + exp(nu[idx] * lchoose(s, k) + k * psi[idx] - lZ) *
        ispos(q[idx] - k + 0.5)
    }
    p[idx] <- ps
  }

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname combinom
#' @export
#' @import RTMB
qcombinom <- function(p, size, prob, nu = 1, lower.tail = TRUE, log.p = FALSE) {

  check_combinom_pars(size, prob, nu)

  np <- max(length(p), length(size), length(prob), length(nu))
  p <- recycle_combinom(p, np)
  size <- recycle_combinom(size, np)
  prob <- recycle_combinom(prob, np)
  nu <- recycle_combinom(nu, np)

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  # tolerate the rounding error of pcombinom(), whose summation is not clamped
  if (any(p < -1e-8 | p > 1 + 1e-8)) stop("p must be in [0,1]")
  p <- pmin(pmax(p, 0), 1)

  psi <- log(prob) - log1p(-prob)

  out <- numeric(np)
  for (i in seq_len(np)) {
    s <- size[i]
    lw <- nu[i] * lchoose(s, 0:s) + (0:s) * psi[i]
    w <- exp(lw - max(lw))
    cdf <- cumsum(w) / sum(w)
    # smallest k in 0:size with F(k) >= p. The fuzz factor (as used by
    # stats::qbinom) makes the q(p(x)) round-trip robust to the rounding
    # difference between this cumulative sum and pcombinom()'s masked sum.
    out[i] <- min(sum(cdf < p[i] * (1 - 64 * .Machine$double.eps)), s)
  }
  out
}

#' @rdname combinom
#' @export
#' @importFrom stats runif
rcombinom <- function(n, size, prob, nu = 1) {

  check_combinom_pars(size, prob, nu)

  size <- recycle_combinom(size, n)
  prob <- recycle_combinom(prob, n)
  nu <- recycle_combinom(nu, n)

  qcombinom(runif(n), size = size, prob = prob, nu = nu)
}
