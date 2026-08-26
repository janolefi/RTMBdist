#' Bell distribution
#'
#' Probability mass function, distribution function, quantile function, and
#' random generation for the Bell distribution.
#'
#' @details
#' This implementation of \code{dbell} and \code{pbell} allows for automatic
#' differentiation with \code{RTMB} with respect to \code{theta}.
#'
#' The Bell distribution (Castellares et al. 2018) is a one-parameter
#' distribution for overdispersed counts with probability mass function
#'
#' \deqn{P(X = x;\, \theta) = \frac{e^{1 - e^{\theta}}\, \theta^{x}\, B_x}{x!}, \quad x = 0, 1, 2, \ldots,}
#'
#' for \eqn{\theta > 0}, where \eqn{B_x} is the \eqn{x}-th Bell number, i.e.
#' the number of ways a set of \eqn{x} elements can be partitioned into
#' non-empty subsets.
#'
#' Its mean and variance are
#' \deqn{E(X) = \theta e^{\theta}, \qquad \mathrm{Var}(X) = \theta e^{\theta} (1 + \theta),}
#' so the dispersion index is \eqn{\mathrm{Var}(X) / E(X) = 1 + \theta > 1}
#' and the distribution is always overdispersed. Note that, unlike the
#' negative binomial or generalised Poisson distributions, the Bell
#' distribution has no separate dispersion parameter: the amount of
#' overdispersion is tied to the mean. As \eqn{\theta \to 0} it approaches the
#' Poisson distribution.
#'
#' The distribution arises as a compound Poisson sum
#' \deqn{X = \sum_{i=1}^{N} Y_i, \qquad N \sim \mathrm{Pois}(e^{\theta} - 1), \quad Y_i \sim \mathrm{ztPois}(\theta),}
#' with the \eqn{Y_i} independent of \eqn{N}, and \code{rbell} uses exactly
#' this representation. It is therefore infinitely divisible, and a member of
#' the one-parameter exponential family with natural parameter \eqn{\log\theta}
#' and sufficient statistic \eqn{x}.
#'
#' The Bell numbers grow faster than any exponential and overflow double
#' precision at \eqn{x = 219}, so \eqn{\log B_x} is used throughout rather than
#' \eqn{B_x}. As \code{x} is data, \eqn{\log B_x} is constant with respect to
#' the parameters and is cached across calls.
#'
#' Neither the distribution function nor the quantile function has a closed
#' form; both are obtained by summing the probability mass function, with
#' \code{qbell} choosing its summation range automatically from the mean and
#' variance.
#'
#' See \code{\link{bell2}} for the parameterisation by the mean.
#'
#' @param x,q integer vector of counts
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param theta vector of positive Bell parameters
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dbell} gives the probability mass function, \code{pbell} gives the distribution function, \code{qbell} gives the quantile function, and \code{rbell} generates random deviates.
#'
#' @references
#' Castellares, F., Ferrari, S. L. P., and Lemonte, A. J. (2018). On the Bell
#' distribution and its associated regression model for count data.
#' \emph{Applied Mathematical Modelling} 56, 172-185.
#' doi:10.1016/j.apm.2017.12.014
#'
#' @examples
#' set.seed(123)
#' x <- rbell(1, 1)
#' d <- dbell(x, 1)
#' p <- pbell(x, 1)
#' q <- qbell(p, 1)
#'
#' # mean and variance
#' theta <- 0.8
#' xs <- 0:100
#' sum(xs * dbell(xs, theta)) # theta * exp(theta)
#' @name bell
NULL

# Recycle to length n, using integer indexing so that advectors are handled.
recycle_bell <- function(x, n) {
  if (length(x) == n) return(x)
  x[rep_len(seq_along(x), n)]
}

# ---------------------------------------------------------------------------
# log Bell numbers
# ---------------------------------------------------------------------------
# B_x appears in the pmf and overflows double precision at x = 219, so only
# log(B_x) is ever formed. Since x is data, this is a constant with respect to
# the parameters: no AD is involved and the values can be cached.

# Largest index served from the cached Bell triangle. The triangle is exact but
# costs O(N^2); beyond this Dobinski's formula is used instead, which is
# accurate to machine precision and cheap for a handful of scattered values.
bell_tri_max <- 20000L

# Cache holding log B_0, ..., log B_N together with the state needed to extend
# it: the current row of the Bell triangle, scaled by exp(-s) to stay in range.
bell_cache <- new.env(parent = emptyenv())
bell_cache$lb <- 0 # lb[n + 1] = log(B_n); initially B_0 = 1 only
bell_cache$row <- 1
bell_cache$s <- 0

# Extend the cache to cover log B_N. Each row of the Bell triangle is the
# cumulative sum of the previous row rotated by one, so a row costs a single
# vectorised cumsum, and B_n is the first entry of row n.
bell_grow <- function(N) {
  cur <- length(bell_cache$lb) - 1L
  if (N <= cur) return(invisible(NULL))

  lb <- c(bell_cache$lb, numeric(N - cur))
  row <- bell_cache$row
  s <- bell_cache$s

  for (n in (cur + 1L):N) {
    m <- length(row)
    row <- cumsum(c(row[m], row))
    lb[n + 1L] <- log(row[1L]) + s
    mx <- row[m + 1L] # the row is increasing, so this is its maximum
    if (mx > 1e250) {
      row <- row / mx
      s <- s + log(mx)
    }
  }

  bell_cache$lb <- lb
  bell_cache$row <- row
  bell_cache$s <- s
  invisible(NULL)
}

# log B_n from Dobinski's formula, log B_n = -1 + log sum_k k^n / k!, summed
# with the log-sum-exp trick. All terms are positive, so no cancellation can
# occur.
#
# The terms t_k = n log k - log k! are sharply peaked: t_k is maximised at the
# saddle point r where r log r = n, that is r = n / W(n), and behaves like a
# Gaussian of width sigma = (n / r^2 + 1 / r)^(-1/2) around it. Summing over a
# window of +/- 40 sigma therefore leaves out only terms of relative size
# exp(-800), while keeping the amount of work per value independent of how far
# the peak has moved out. Summing from k = 1 instead would cost O(n).
# Vectorised over n. Consecutive windows overlap almost entirely, so log(k)
# and log(k!) are formed once over the union of all the windows and then only
# sliced, which is what keeps the cost per value flat.
lbell_dobinski <- function(n) {
  out <- numeric(length(n))
  small <- n <= 1
  if (all(small)) return(out)

  ns <- n[!small]
  r <- ns / lambertW.num(ns)
  half <- pmax(60, ceiling(40 / sqrt(ns / r^2 + 1 / r)))
  lo <- pmax(1, floor(r - half))
  hi <- ceiling(r + half)

  gl <- min(lo)
  k <- gl:max(hi)
  lk <- log(k)
  lf <- lgamma(k + 1)

  res <- numeric(length(ns))
  for (i in seq_along(ns)) {
    j <- (lo[i] - gl + 1):(hi[i] - gl + 1)
    lt <- ns[i] * lk[j] - lf[j]
    m <- max(lt)
    res[i] <- m + log(sum(exp(lt - m))) - 1
  }
  out[!small] <- res
  out
}

# log B_n for a vector of non-negative integers.
lbell <- function(n) {
  n <- as.numeric(n)
  out <- numeric(length(n))
  big <- n > bell_tri_max
  if (any(!big)) {
    bell_grow(max(n[!big]))
    out[!big] <- bell_cache$lb[n[!big] + 1]
  }
  if (any(big)) out[big] <- lbell_dobinski(n[big])
  out
}

# ---------------------------------------------------------------------------

# Shared argument check.
check_bell_pars <- function(theta) {
  if (any(theta <= 0)) stop("theta must be > 0")
  if (any(!is.finite(theta))) stop("theta must be finite")
  invisible(NULL)
}

# log-pmf on the grid 0:N for a single theta, which may be an advector. The
# k = 0 entry is written out separately so that no 0 * log(theta) is formed,
# which would be NaN rather than 0 should theta ever underflow to zero.
lpmf_bell <- function(N, theta) {
  k <- 0:N
  out <- k * log(theta) - expm1(theta) + lbell(k) - lgamma(k + 1)
  out[1] <- -expm1(theta)
  out
}

#' @rdname bell
#' @export
#' @import RTMB
dbell <- function(x, theta, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    check_bell_pars(theta)
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dbell", x = x, theta = theta, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dbell", x = x, theta = theta, log = log))
  }
  if (inherits(x, "advector")) {
    stop("dbell is not differentiable with respect to x, because the Bell numbers are only defined on the integers. Use method = \"cdf\" for OSA residuals.")
  }

  nx <- max(length(x), length(theta))
  x <- recycle_bell(x, nx)
  theta <- recycle_bell(theta, nx)

  xn <- as.numeric(x)
  ok <- is.finite(xn) & xn >= 0 & xn == floor(xn)
  xs <- ifelse(ok, xn, 0) # placeholder; invalid entries become -Inf below

  logdens <- xs * log(theta) - expm1(theta) + lbell(xs) - lgamma(xs + 1)

  # x = 0 contributes no log(theta) term at all; assigning it directly avoids
  # 0 * log(theta), which is NaN rather than 0 if theta underflows to zero.
  z <- xs == 0
  logdens[z] <- -expm1(theta[z])
  logdens[!ok] <- -Inf # non-integer and negative x carry no mass

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname bell
#' @export
#' @import RTMB
pbell <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    check_bell_pars(theta)
  }

  if (inherits(q, "advector")) {
    stop("pbell is not differentiable with respect to q, because the Bell numbers are only defined on the integers.")
  }

  nq <- max(length(q), length(theta))
  q <- recycle_bell(q, nq)
  theta <- recycle_bell(theta, nq)

  qq <- floor(as.numeric(q))
  inrange <- which(qq >= 0) # which() drops the NAs, which are left at p = 0

  if (inherits(theta, "advector")) {
    # One accumulation per observation. theta varies by observation in an AD
    # context, so there is nothing to share between them.
    p <- 0 * theta
    for (i in inrange) p[i] <- sum(exp(lpmf_bell(qq[i], theta[i])))
  } else {
    # theta is data, so all entries sharing a value share one cumulative sum.
    p <- numeric(nq)
    for (th in unique(theta[inrange])) {
      idx <- inrange[theta[inrange] == th]
      cdf <- cumsum(exp(lpmf_bell(max(qq[idx]), th)))
      p[idx] <- cdf[qq[idx] + 1]
    }
    p <- pmin(pmax(p, 0), 1) # a CDF, up to rounding in the summation
  }

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname bell
#' @export
qbell <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {

  check_bell_pars(theta)

  np <- max(length(p), length(theta))
  p <- rep_len(as.numeric(p), np)
  theta <- rep_len(as.numeric(theta), np)

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (any(!is.na(p) & (p < 0 | p > 1))) stop("p must be in [0,1]")

  out <- numeric(np)
  out[is.na(p)] <- NA_real_
  # as for qpois(), the upper endpoint of the support is unbounded
  out[which(p >= 1)] <- Inf
  todo <- which(p < 1) # which() drops the NAs

  for (th in unique(theta[todo])) {
    idx <- todo[theta[todo] == th]
    ptop <- max(p[idx])

    # start from the mean plus ten standard deviations and extend if the
    # requested probability lies further out in the tail
    mu <- th * exp(th)
    N <- max(10, ceiling(mu + 10 * sqrt(mu * (1 + th))))
    repeat {
      cdf <- cumsum(exp(lpmf_bell(N, th)))
      if (cdf[N + 1] >= ptop || cdf[N + 1] > 1 - 1e-14 || N >= 1e7) break
      N <- 2 * N
    }

    # smallest k with F(k) >= p. The fuzz factor, as used by stats::qpois(),
    # keeps the q(p(x)) round-trip robust to rounding in the summation.
    out[idx] <- vapply(p[idx],
                       function(pp) sum(cdf < pp * (1 - 64 * .Machine$double.eps)),
                       numeric(1))
  }

  out
}

#' @rdname bell
#' @export
#' @importFrom stats rpois
rbell <- function(n, theta) {

  check_bell_pars(theta)

  if (length(n) > 1) n <- length(n)
  theta <- rep_len(as.numeric(theta), n)

  # compound Poisson representation: a Poisson number of zero-truncated
  # Poisson variates. Exact, with no truncation of the support.
  N <- rpois(n, expm1(theta))

  out <- numeric(n)
  tot <- sum(N)
  if (tot > 0) {
    grp <- rep.int(seq_len(n), N)
    y <- rztpois(tot, theta[grp])
    s <- rowsum(y, group = grp, reorder = FALSE)
    out[as.integer(rownames(s))] <- as.numeric(s)
  }
  out
}
