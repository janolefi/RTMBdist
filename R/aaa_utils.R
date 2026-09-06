# getting OSA residual and simulation functions from RTMB (not exported)
dGenericOSA <- get("dGenericOSA", envir = asNamespace("RTMB"), inherits = FALSE)
dGenericSim <- get("dGenericSim", envir = asNamespace("RTMB"), inherits = FALSE)

# getting ad_context from RTMB (not exported)
ad_context <- get("ad_context", envir = asNamespace("RTMB"), inherits = FALSE)

#' AD-compatible error function and complementary error function
#'
#' @param x vector of evaluation points
#'
#' @returns \code{erf(x)} returns the error function and \code{erfc(x)} returns the complementary error function.
#'
#' @examples
#' erf(1)
#' erfc(1)
#' @name erf
NULL
#' @rdname erf
#' @export
#' @importFrom RTMB pnorm
erf <- function(x) {
  2 * RTMB::pnorm(x * sqrt(2)) - 1 # + eps
}
#' @rdname erf
#' @export
erfc <- function(x) {
  1 - erf(x) # + eps
}

#' Lambert W function (principal branch)
#'
#' Solves \eqn{W(x) e^{W(x)} = x} for the principal branch \eqn{W_0}.
#'
#' @param x vector of evaluation points, \eqn{x \ge -1/e}.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' The value is obtained by Halley iteration, which converges to machine
#' precision in a handful of steps over the whole domain. For AD, the function
#' is registered as an atomic operation via \code{\link[RTMB]{ADjoint}} with
#' the analytic derivative
#' \deqn{W'(x) = \frac{1}{e^{W(x)} (1 + W(x))},}
#' expressed through the returned value rather than through \eqn{x}. Written
#' this way the derivative is itself an AD-able expression, so derivatives of
#' every order are available; in particular the third-order derivatives that
#' the gradient of a Laplace approximation requires.
#'
#' The principal branch is defined for \eqn{x \ge -1/e}, with
#' \eqn{W_0(-1/e) = -1} and \eqn{W_0(x) \ge -1} throughout. Values below
#' \eqn{-1/e} return \code{NaN} with a warning. The derivative is infinite at
#' the branch point \eqn{x = -1/e}.
#'
#' @returns The principal branch of the Lambert W function evaluated at \code{x}.
#' @export
#'
#' @examples
#' lambertW(exp(1)) # 1
#' lambertW(0) # 0
#' x <- c(0.5, 1, 10, 1000)
#' lambertW(x) * exp(lambertW(x)) - x # ~ 0
lambertW <- function(x) {
  if (inherits(x, "advector")) return(lambertW.ad(x))
  lambertW.num(x)
}

# Double-precision principal branch, by Halley iteration. Vectorised, and free
# of NaN warnings: each initial-guess branch is only ever evaluated on the
# subset of x it is valid for.
lambertW.num <- function(x) {
  x <- as.numeric(x)
  if (any(x < -exp(-1), na.rm = TRUE))
    warning("lambertW: x < -1/e is outside the principal branch")

  y <- rep(NaN, length(x))
  ok <- !is.na(x) & x >= -exp(-1)
  xo <- x[ok]

  # Initial guess, by regime. Each branch is only ever evaluated on the subset
  # it is valid for, so no NaNs are produced along the way.
  yo <- xo / (1 + xo)
  near <- xo < -0.3 # branch point: series in p = sqrt(2 (e x + 1))
  p <- sqrt(2 * pmax(exp(1) * xo[near] + 1, 0))
  yo[near] <- -1 + p - p^2 / 3 + 11 * p^3 / 72
  hi <- xo > 1 # asymptotic W(x) ~ log(x) - log(log(x))
  yo[hi] <- log(xo[hi]) - log(log1p(xo[hi]))

  # Halley. Converges in ~3 steps; 12 leaves a wide margin near the branch
  # point, where convergence is slowest and the denominator vanishes.
  for (i in 1:12) {
    e <- exp(yo)
    f <- yo * e - xo
    step <- f / (e * (yo + 1) - (yo + 2) * f / (2 * yo + 2))
    step[!is.finite(step)] <- 0
    yo <- pmax(yo - step, -1) # the principal branch never goes below -1
  }
  y[ok] <- yo
  y
}

# AD version. df is written in terms of the returned value y (and uses only
# AD-able operations), so RTMB can differentiate it again, and again: the
# recursion closes and derivatives of all orders are exact. Using
# 1 / (exp(y) (1 + y)) rather than the equivalent y / (x (1 + y)) also keeps
# the derivative finite at x = 0, where the latter is 0/0.
lambertW.ad <- RTMB::ADjoint(
  f = function(x) lambertW.num(x),
  df = function(x, y, dy) dy / (exp(y) * (1 + y)),
  name = "lambertW"
)

#' Smooth approximation to the absolute value function
#'
#' @param x vector of evaluation points
#' @param epsilon smoothing constant
#'
#' @details
#' We approximate the absolute value here as
#' \deqn{\vert x \vert \approx \sqrt{x^2 + \epsilon}}
#'
#' @returns Smooth absolute value of \code{x}.
#' @export
#'
#' @examples
#' abs(0)
#' abs_smooth(0, 1e-4)
abs_smooth <- function(x, epsilon = 1e-6) {
  sqrt(x^2 + epsilon)
}

## AD pmin/pmax helpers that work for both ad and numeric:
pmin.ad <- function(x, y) apply(cbind(x,y), 1, min)
pmax.ad <- function(x, y) apply(cbind(x,y), 1, max)

# Replace with loop versions for now to avoid method dispatch error in r-devel
# pmin.ad <- function(x, y) {
#   n <- max(length(x), length(y))
#   x <- x[rep_len(seq_along(x), n)]   # recycle via integer indexing;
#   y <- y[rep_len(seq_along(y), n)]   # "[" has an explicit advector method
#   out <- x
#   for (j in seq_len(n)) out[j] <- min(x[j], y[j])
#   out
# }
# pmax.ad <- function(x, y) {
#   n <- max(length(x), length(y))
#   x <- x[rep_len(seq_along(x), n)]
#   y <- y[rep_len(seq_along(y), n)]
#   out <- x
#   for (j in seq_len(n)) out[j] <- min(-x[j], -y[j]) * (-1)
#   out
# }

## AD-indicator constructors
# 1 if x == 0, 0 otherwise
iszero <- function(x) {
  if(inherits(x, c("advector", "osa", "simref"))) {
    return(iszero.ad(x))
  } else {
    return(as.numeric(x == 0))
  }
}

iszero.ad <- RTMB::ADjoint(f = function(x) as.numeric(x==0),
                           df = function(x,y,dy) RTMB::AD(rep(0, length(x))),
                           name = "iszero.ad")
# zero <- ADjoint(f = function(x) rep(0, length(x)),
#                 df = function(x, y, dy) zero(x),
#                 name = "zero")
# iszero <- ADjoint(f = function(x) as.numeric(x==0),
#                   df = function(x,y,dy) zero(x),
#                   name = "iszero")
# 1 if x != 0, 0 otherwise
isnonzero <- function(x) {
  1 - iszero(x)
}
# 1 if x => 0, 0 otherwise
ispos <- function(x) {
  s <- sign(x)
  0.5 * (s + abs(s))
}
# 1 if x < 0, 0 otherwise
isneg <- function(x) {
  s <- sign(x)
  - 0.5 * (s - abs(s))
}
# 1 if x > 0, 0 otherwise
ispos_strict <- function(x) {
  ispos(x) * isnonzero(x)
}
# 1 if x < val, 0 otherwise
smaller <- function(x, val) {
  s <- sign(x - val)
  0.5 * (abs(s) - s)
}
# 1 if x > val, 0 otherwise
greater <- function(x, val) {
  s <- sign(val - x)
  0.5 * (abs(s) - s)
}
# turns +/-Inf into largest finite value
as.finite <- function(x) {
  x <- pmin.ad(x, .Machine$double.xmax)
  x <- pmax.ad(x, -.Machine$double.xmax)
  return(x)
}
as.finite.neg <- function(x) {
  pmax.ad(x, -.Machine$double.xmax)
}


## Logarithm of zero-inflated density/ pmf
# x == 0: p0
# x > 0: (1-p0) * pdf(x)
log_zi <- function(x, logdens, zeroprob) {
  logdens <- as.finite(logdens) # turn +/- Inf into finite
  logdens <- RTMB::logspace_add(
    log(iszero(x)) + log(zeroprob),
    log(isnonzero(x)) + log1p(-zeroprob) + logdens
  )
}
# x == 0: p0 + pmf(0)
# x > 0: (1-p0) * pmf(x)
log_zi_discrete <- function(x, logdens, zeroprob) {
  RTMB::logspace_add(
    log(iszero(x)) + log(zeroprob),
    log1p(-zeroprob) + logdens
    )
}
# log Beta function
lbeta.ad <- function(a, b) {
  # lgamma(a) + lgamma(b) - lgamma(a + b)
  lbeta(a,b)
}

# Log multivariate gamma, AD-friendly
lmultigamma <- function(a, p) {
  # Only check bounds if not in AD context
  if (!ad_context()) {
    if (a <= (p - 1) / 2) stop("a must be greater than (p - 1) / 2")
  }
  sum(lgamma(a + (1 - 1:p)/2))
}

# Inverse Box-Cox transformation shared by the BCCG, BCT and BCPE quantile
# functions: maps z on the standardised scale back to the response scale.
# At nu = 0 the first branch evaluates to 1^Inf = 1, so ifelse() stays finite.
inv_boxcox <- function(mu, sigma, nu, z) {
  ifelse(nu != 0, mu * (nu * sigma * z + 1)^(1 / nu), mu * exp(sigma * z))
}

reggamma <- function(s, x) {
  if (ad_context()) {
    # RTMB's pgamma errors when lower.tail is passed inside an AD context, so
    # the upper tail is formed by complement there; outside AD the direct
    # upper-tail evaluation is kept because it is accurate far into the tail
    return(1 - pgamma(x, shape = s, scale = 1))
  }

  pgamma(x, shape = s, scale = 1, lower.tail = FALSE)
}

# Generate helpful error message if user wrote likeliood in wrong order to simulate
simulation_check <- function(args, exclude = c("x", "log")) {
  args <- args[setdiff(names(args), exclude)]
  if (any(vapply(args, function(a) inherits(a, "simref"), logical(1)))) {
    stop(
      "Automatic simulation requires the likelihood to follow the model hierarchy: random effects first, then data given those random effects."
    )
  }
}
