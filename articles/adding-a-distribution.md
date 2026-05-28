# Guide to adding a distribution

## Background

RTMBdist provides additional probability distributions that work inside
[`RTMB::MakeTape`](https://rdrr.io/pkg/RTMB/man/Tape.html) for automatic
differentiation (AD). RTMB records a computation tape by operator
overloading: every arithmetic operation on an `advector` is intercepted.
All code must therefore be plain R arithmetic — anything that drops
below R to compiled C++ is invisible to the tape. This is why built-in
distributions such as `dgamma` must be called as `RTMB::dgamma`, not
[`stats::dgamma`](https://rdrr.io/r/stats/GammaDist.html), and why
`@import RTMB` is needed whenever a function uses mathematical
operations (`exp`, `log`, `sqrt`, etc.).

Every density function must support three modes:

1.  **Plain R evaluation** — standard use outside of RTMB.
2.  **AD tape recording** — called with an `advector`; gradients must
    flow.
3.  **Simulation and OSA residuals** — called with a `simref` or `osa`
    object.

## What to implement

At minimum, provide a **density** (`d<dist>`) and a **random-number
generator** (`r<dist>`). A **CDF** (`p<dist>`) and **quantile function**
(`q<dist>`) are strongly preferred; they are required for OSA residuals
and quantile-based tests. All four functions must be pure R arithmetic.

## Naming conventions

Follow the `d`/`p`/`q`/`r` prefix convention. Spell out distribution
names rather than abbreviating — `gompertz`, `invgamma`, `laplace` —
unless the name becomes unwieldy. Reparameterised variants take a `2`
suffix: `beta2`, `gamma2`.

## Files to create or modify

| File                           | What to do                             |
|--------------------------------|----------------------------------------|
| `R/<dist>.R`                   | implement `d`, `p`, `q`, `r` functions |
| `tests/testthat/test-<dist>.R` | add standard tests                     |
| `vignettes/distlist.Rmd`       | add one bullet in alphabetical order   |

------------------------------------------------------------------------

## Implementing the density function

``` r

dfoo <- function(x, theta1, theta2 = 1, log = FALSE) {

  # 1. Validate inputs — skip during AD recording
  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args)           # catches likelihood written in wrong order
    if (any(theta1 <= 0)) stop("theta1 must be positive")
  }

  # 2. Simulation and OSA escapes
  if (inherits(x, "simref"))
    return(dGenericSim("dfoo", x = x, theta1 = theta1, theta2 = theta2, log = log))
  if (inherits(x, "osa"))
    return(dGenericOSA("dfoo", x = x, theta1 = theta1, theta2 = theta2, log = log))

  # 3. Log-density
  logdens <- ...

  if (log) return(logdens)
  exp(logdens)
}
```

A few rules:

- Always compute in log space and only call
  [`exp()`](https://rdrr.io/r/base/Log.html) at the very end.
- The string passed to `dGenericSim`/`dGenericOSA` must exactly match
  the function name, because RTMB uses it to look up the corresponding
  `p<dist>` or `r<dist>`.
- Parameter validation and `simulation_check` must be inside the
  `if (!ad_context())` guard — they would fail or give wrong results on
  `advector` inputs.

## Implementing the CDF

``` r

pfoo <- function(q, theta1, theta2 = 1, lower.tail = TRUE, log.p = FALSE) {
  if (!ad_context()) {
    if (any(theta1 <= 0)) stop("theta1 must be positive")
  }

  p <- ...   # compute CDF using RTMB:: functions where needed

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}
```

No `simulation_check` is needed in `p<dist>`. Apply `lower.tail` and
`log.p` at the end, in that order.

## Implementing the quantile function

``` r

qfoo <- function(p, theta1, theta2 = 1, lower.tail = TRUE, log.p = FALSE) {
  if (!ad_context()) {
    if (any(theta1 <= 0)) stop("theta1 must be positive")
  }

  if (log.p) p <- exp(p)          # undo transformations first ...
  if (!lower.tail) p <- 1 - p     # ... in the reverse order of pfoo

  ...   # compute quantile
}
```

## Implementing the RNG

``` r

rfoo <- function(n, theta1, theta2 = 1) {
  if (any(theta1 <= 0)) stop("theta1 must be positive")
  ...
}
```

The RNG is never called during AD, so no `ad_context()` guard is needed.
Two common approaches:

- **Inversion**: `qfoo(runif(n), theta1, theta2)` — works whenever
  `q<dist>` exists.
- **Direct formula**: use a known relationship to a base distribution,
  e.g. `log(1 + rexp(n, rate = eta)) / b` for the Gompertz distribution.

------------------------------------------------------------------------

## Numerical stability

Prefer numerically stable forms whenever the argument may be close to
zero:

| Instead of                   | Write                |
|------------------------------|----------------------|
| `exp(x) - 1`                 | `expm1(x)`           |
| `log(1 + x)`                 | `log1p(x)`           |
| `1 - exp(-x)`                | `-expm1(-x)`         |
| `x^2`, `x^3`                 | `x * x`, `x * x * x` |
| `x^a` (non-integer exponent) | `exp(a * log(x))`    |

Avoid expressions that produce `0 * (-Inf)` — this evaluates to `NaN`
and kills the gradient. Common fixes:

- Add a small offset before taking a log:
  `log(x + .Machine$double.xmin)`.
- Clamp a log term that can hit `-Inf` at a boundary:
  `as.finite(log(...))`.

## Conditional logic

Never branch on the value of `x` inside the tape with plain `if/else` or
`ifelse` — the condition is not recorded and the gradient will be wrong.
Use the smooth indicator helpers from `aaa_utils.R` instead:

| Helper         | Meaning              |
|----------------|----------------------|
| `iszero(x)`    | 1 if x == 0, else 0  |
| `isnonzero(x)` | 1 if x != 0, else 0  |
| `ispos(x)`     | 1 if x \>= 0, else 0 |
| `isneg(x)`     | 1 if x \< 0, else 0  |

These are defined as pure arithmetic and propagate gradients correctly.
Look distributions already implemented to see how to use this.

------------------------------------------------------------------------

## Zero-inflated distributions

Two utility functions handle the mixture in a tape-compatible way:

- `log_zi(x, logdens, zeroprob)` — continuous ZI: point mass at 0,
  density for x \> 0.
- `log_zi_discrete(x, logdens, zeroprob)` — discrete ZI: adds inflation
  mass on top of the PMF at 0.

Typical usage:

``` r

logdens <- RTMB::dgamma(x + .Machine$double.xmin, shape = shape, scale = scale, log = TRUE)
logdens <- log_zi(x, logdens, zeroprob)
```

The `eps` offset prevents `log(0)` before `log_zi` has a chance to zero
out the weight.

## Reparameterised variants

Transform parameters at the top and delegate to the base function:

``` r

dfoo2 <- function(x, mean, sd, log = FALSE) {
  if (!ad_context()) { ... }
  if (inherits(x, "simref")) return(dGenericSim("dfoo2", ...))
  if (inherits(x, "osa"))    return(dGenericOSA("dfoo2", ...))

  shape <- mean * mean / (sd * sd)
  scale <- sd * sd / mean
  dfoo(x, shape = shape, scale = scale, log = log)
}
```

Do not skip the `simulation_check` and simref/osa escapes just because
the base function also has them.

------------------------------------------------------------------------

## Documentation

Use the `@name` / `NULL` trick to attach the shared Roxygen block to a
dummy object and give each function `@rdname`:

``` r

#' Foo distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Foo distribution.
#'
#' @details
#' The Foo distribution with parameter \eqn{\theta > 0} has density
#' \deqn{f(x;\,\theta) = \theta e^{-\theta x}, \quad x \ge 0.}
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param theta positive rate parameter
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities are
#'   returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @references \url{https://en.wikipedia.org/wiki/...}
#' @return \code{dfoo} gives the density, \code{pfoo} gives the distribution
#'   function, \code{qfoo} gives the quantile function, and \code{rfoo}
#'   generates random deviates.
#' @examples
#' x <- rfoo(1, theta = 2)
#' d <- dfoo(x, theta = 2)
#' p <- pfoo(x, theta = 2)
#' q <- qfoo(p, theta = 2)
#' @name foo
NULL

#' @rdname foo
#' @export
#' @import RTMB
dfoo <- function(x, theta, log = FALSE) { ... }
```

Write the `\deqn{}` in clean mathematical form — not the numerically
stable implementation. Use `@import RTMB` on every function that uses
arithmetic or mathematical operations. Use `@importFrom stats r<base>`
in the RNG for the base sampler.

------------------------------------------------------------------------

## Tests

Place tests in `tests/testthat/test-<dist>.R`. Use two parameter
combinations and always include an AD gradient check.

**Continuous distribution:**

``` r

test_that("foo passes standard distribution checks (theta=1)", {
  check_continuous_dist(dfoo, pfoo, qfoo,
    xs = c(0.2, 0.5, 1, 2),
    lower = 0, upper = Inf,
    theta = 1)
})

test_that("foo passes standard distribution checks (theta=3)", {
  check_continuous_dist(dfoo, pfoo, qfoo,
    xs = c(0.1, 0.4, 0.8, 1.5),
    lower = 0, upper = Inf,
    theta = 3)
})

test_that("foo AD gradient has no NaN", {
  check_ad_gradient(dfoo, rfoo, theta = 1)
})
```

Choose `xs` so that `pfoo(xs)` stays well away from 0 and 1 — the
`q(p(x))` round-trip returns `Inf` or `-Inf` when the CDF saturates in
double precision.

Other helpers for different distribution types:

| Type | Helper |
|----|----|
| Continuous | `check_continuous_dist` |
| Continuous without quantile function | `check_continuous_dist(..., qfun = NULL)` |
| Zero-inflated continuous | `check_zeroinfl_dist` |
| Mixed/inflated (ZI+OI etc.) | `check_inflated_dist` |
| Discrete | `check_discrete_dist` |

## Vignette entry

Add one bullet to `vignettes/distlist.Rmd` in the correct section
(Continuous / Discrete / Multivariate) in alphabetical order:

    * [`foo(theta)`](../reference/foo.html): Foo distribution parameterised by rate `theta`
