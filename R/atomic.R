# Atomic tapes: a function taped once with RTMB::MakeTape() and turned into a single
# operation with $atomic(). Each call then adds one node to the calling tape instead of
# all operations of the function, while derivatives of all orders still work.
# The tapes are external pointers, so they cannot be stored in the namespace when the package
# is built; they are created on first use and cached for the session.
.atomic_cache <- new.env(parent = emptyenv())

# the atomic tape of the scalar function f(p), created at the point x on first use
get_atomic <- function(name, f, x) {
  A <- .atomic_cache[[name]]
  if (is.null(A)) {
    A <- RTMB::MakeTape(f, x)$atomic()
    assign(name, A, envir = .atomic_cache)
  }
  A
}

# apply the atomic A elementwise to the recycled arguments in ..., one call per element
map_atomic <- function(A, ...) {
  args <- list(...)
  n <- max(lengths(args))
  args <- lapply(args, function(x) x[(seq_len(n) - 1) %% length(x) + 1])
  # AD() on each input: within the package namespace, c() of a plain number followed by an
  # AD variable drops the AD class (an advector is stored as a complex vector)
  do.call(c, lapply(seq_len(n), function(i) A(do.call(c, lapply(args, function(x) AD(x[i]))))))
}

# Regularised lower incomplete gamma function P(shape, x) = pgamma(x, shape), for use in
# distribution functions. Outside AD this is stats::pgamma(). In AD context, RTMB's pgamma()
# is avoided for small x, where its second and third derivatives in x are not finite
# (x < 1, and NaN in the second derivative at x = 1 exactly; RTMB 2.0). These are needed by
# sdreport() and the Laplace approximation. For x < 2 the series
# P(a, x) = x^a e^(-x) sum_j x^j / ((a + 1) ... (a + j)) / Gamma(a + 1) is used instead; its
# terms are positive, and 24 of them give a relative truncation error below 1e-16 for any
# shape. For x <= 0 the result is 0.
pgamma_ad <- function(x, shape) {
  if (!ad_context()) return(stats::pgamma(x, shape = shape))
  map_atomic(get_atomic("pgamma_ad", pgamma_ad_body, c(1, 1)), x, shape)
}

pgamma_ad_body <- function(p) {
  x <- p[1]
  a <- p[2]
  x0 <- 2
  pos <- 1 - smaller(x, 1e-300) # 1 if x > 0
  near <- pos * smaller(x, x0)
  far <- pos - near

  # Each branch is evaluated at a harmless point where it is not used, so that no infinite
  # or NaN derivative can reach the result (0 * NaN would still be NaN): the series at x = 1,
  # which is fine for it, and RTMB's pgamma() beyond x0. as.finite() keeps 0 * x finite for
  # x = +-Inf.
  xfin <- as.finite(x)
  xn <- near * xfin + (1 - near)
  term <- 1
  M <- 1
  for (j in 1:24) {
    term <- term * xn / (a + j)
    M <- M + term
  }
  series <- exp(a * log(xn) - xn - lgamma(a + 1)) * M

  xf <- far * xfin + (1 - far) * 2 * x0
  near * series + far * RTMB::pgamma(xf, shape = a, scale = 1)
}
