# Distribution function of a discrete distribution on 0, 1, 2, ... by summing its probability
# mass function, AD-compatible in the parameters. Used by p-functions of discrete
# distributions whose distribution function has no closed form.
#
# dfun must be vectorised and AD-compatible in its parameters, which are passed as a named
# list `pars` (a list rather than ..., so that no parameter name can partially match an
# argument) and recycled together with q. q must be numeric: it decides how many terms are
# summed, which could not change when a tape is re-evaluated. This is how RTMB's CDF-based
# OSA residuals use it, as they evaluate the distribution function at the data.
# For a bounded support 0, ..., upper, `upper` must be numeric as well. Both tails are then
# summed directly and divided by the sum over the whole support, so that p is exactly 1 at
# q >= upper (rounding could otherwise put it above 1, and log(1 - p) in the OSA residuals
# would be NaN) and small upper tails stay accurate. For an unbounded support, p is capped
# to [0, 1] for the same reason, and the cost grows with max(q).
discrete_cdf <- function(dfun, q, pars, upper = Inf, lower.tail = TRUE, log.p = FALSE) {
  if (inherits(q, "advector") || inherits(upper, "advector")) {
    stop("q and the end of the support must be numeric data, not AD variables.")
  }

  # names of the parameters would otherwise end up on the result
  pars <- lapply(pars, function(x) if (inherits(x, "advector")) x else unname(x))
  n <- max(length(q), lengths(pars), length(upper))
  shared <- all(lengths(pars) == 1) && length(unique(upper)) == 1
  rec <- function(x) x[(seq_len(n) - 1) %% length(x) + 1]
  q <- floor(rec(q))
  upper <- rec(upper)
  bounded <- is.finite(upper)
  if (any(bounded) && !all(bounded)) stop("upper must be finite for all elements or for none.")
  qc <- pmin(q, upper)

  # the sum runs up to the largest q, or over the whole support where it is bounded
  last <- ifelse(bounded, upper, qc)
  last <- last[is.finite(last)]
  K <- if (length(last)) max(0, last) else 0

  if (shared) {
    # one set of parameters: the mass function is evaluated once and summed cumulatively
    d <- do.call(dfun, c(list(0:K), pars))
    zero <- 0 * d[1]
    below <- c(zero, cumsum(d))       # below[j + 2] = P(X <= j), below[1] = 0
    above <- c(rev(cumsum(rev(d))), zero) # above[j + 2] = P(X > j)
    i <- pmin(pmax(qc, -1), K) + 2
    i[is.na(i)] <- 1
    total <- below[K + 2]
    below <- below[i]
    above <- above[i]
  } else {
    # parameters per element: loop over the support, vectorised over the elements
    pars <- lapply(pars, rec)
    below <- above <- total <- 0
    for (k in 0:K) {
      # pmin() keeps dfun inside a bounded support, where the term is then masked out
      pk <- do.call(dfun, c(list(pmin(k, upper)), pars)) * (k <= upper)
      below <- below + pk * (k <= qc & !is.na(qc))
      above <- above + pk * (k > qc & !is.na(qc))
      total <- total + pk
    }
  }

  if (any(bounded)) {
    denom <- total * bounded + (1 - bounded)
    p <- if (lower.tail) below / denom else above / denom
  } else {
    p <- if (lower.tail) below else 1 - below
  }
  # unbounded support: keep p within [0, 1] against rounding, and q = Inf is certain
  p <- p - greater(p, 1) * (p - 1)
  p <- p * (1 - smaller(p, 0))
  inf <- which(!bounded & q == Inf)
  if (length(inf)) p[inf] <- if (lower.tail) 1 else 0

  if (log.p) p <- log(p)
  if (anyNA(q)) p[is.na(q)] <- NA
  p
}
