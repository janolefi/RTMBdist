# Numerical distribution function by integrating a density with RTMB's integrate(),
# which is AD-compatible from RTMB 2.0 on. Used by p-functions in AD context whose
# closed-form or package CDF is not AD-compatible.
#
# dfun must be AD-compatible, vectorised in its first argument, have a `log` argument and
# be positive at q. Its parameters are passed as a named list `pars` (a list rather than
# ..., so that no parameter name can partially match an argument, e.g. df and dfun) and
# recycled together with q, so each element gets its own integral.
# To avoid cancellation, the smaller tail is integrated directly: from `lower` to q if
# q <= centre, and from q to `upper` otherwise, taking the complement if needed.
# The tail is integrated relative to the density at q, i.e. as f(q) * int f(x) / f(q) dx
# with the ratio computed on the log scale. The integral is then of order one also far out
# in the tails, which avoids underflow (subnormal values give NaN derivatives) and keeps
# log.p accurate where the probability itself underflows.
# `centre` should be close to the mode, so that the peak of the density sits at an end of
# the integration range, where the adaptive quadrature resolves it reliably. `scale` sets
# the unit in which infinite tails are integrated and should be of the order of the spread
# of the density. Only the values of `centre` and `scale` at tape time are used. The
# choices they make therefore stay fixed when the tape is re-evaluated at other parameters,
# which leaves the result correct and only affects numerical accuracy.
numerical_cdf <- function(dfun, q, pars, centre, scale = 1, lower = -Inf, upper = Inf,
                          lower.tail = TRUE, log.p = FALSE, rel.tol = 1e-8) {
  n <- max(length(q), lengths(pars))
  rec <- function(x) x[(seq_len(n) - 1) %% length(x) + 1]
  q <- rec(q)
  pars <- lapply(pars, rec)

  qv <- value_of(q)
  left <- qv <= rec(value_of(centre))
  left[is.na(left)] <- TRUE # NA in q, filled in at the end
  h <- rec(value_of(scale))

  # Log of the integral over one tail for scalar q, scale and parameters. RTMB's Vectorize()
  # tapes this once and maps the tape over all elements, which makes tape construction much
  # faster. This also fixes all branching to that of the first element, hence one call per tail.
  # RTMB's AD integrate() can give a zero derivative with respect to a finite integration
  # limit, so the limits are kept fixed and q is moved into the integrand instead.
  arg_names <- c(".q", ".h", names(pars))
  log_tail <- function(from_lower) {
    fun <- function() {
      args <- mget(arg_names)
      q <- args[[1]]
      h <- args[[2]]
      logf <- function(x) do.call(dfun, c(list(x), args[-(1:2)], list(log = TRUE)))
      logf_q <- logf(q)
      r <- function(x) exp(logf(x) - logf_q)
      # abs.tol = 0 so that the relative tolerance decides
      I <- function(g, b) integrate(g, 0, b, rel.tol = rel.tol, abs.tol = 0)$value
      logf_q + log(
        if (from_lower) {
          if (is.infinite(lower)) h * I(function(v) r(q - h * v), Inf)
          else (q - lower) * I(function(u) r(lower + u * (q - lower)), 1)
        } else {
          if (is.infinite(upper)) h * I(function(v) r(q + h * v), Inf)
          else (upper - q) * I(function(u) r(q + u * (upper - q)), 1)
        }
      )
    }
    formals(fun) <- stats::setNames(rep(alist(x = ), length(arg_names)), arg_names)
    Vectorize(fun)
  }

  # infinite q is left out of the integration (it would be a taped variable inside
  # Vectorize), its directly integrated tail is empty
  ls <- rep(-Inf, n)
  if (ad_context()) ls <- advector(ls)
  for (from_lower in c(TRUE, FALSE)) {
    i <- which(left == from_lower & is.finite(qv))
    if (length(i)) {
      ls[i] <- do.call(log_tail(from_lower), c(list(q[i], h[i]), lapply(pars, `[`, i)))
    }
  }

  flip <- left != lower.tail
  p <- if (log.p) ls else exp(ls)
  if (any(flip)) p[flip] <- if (log.p) log1p(-exp(ls[flip])) else 1 - exp(ls[flip])
  if (anyNA(qv)) p[is.na(qv)] <- NA
  p
}
