# Beta-binomial distribution

Density, distribution function and random generation for the
beta-binomial distribution.

## Usage

``` r
dbetabinom(x, size, shape1, shape2, log = FALSE)

pbetabinom(q, size, shape1, shape2, lower.tail = TRUE, log.p = FALSE)

rbetabinom(n, size, shape1, shape2)
```

## Arguments

- x:

  vector of non-negative counts.

- size:

  vector of total counts (number of trials). Needs to be \>= `x`.

- shape1:

  positive shape parameter 1 of the Beta prior.

- shape2:

  positive shape parameter 2 of the Beta prior.

- log:

  logical; if `TRUE`, densities are returned on the log scale.

- q:

  vector of quantiles.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le q\]\\,
  otherwise \\P\[X \> q\]\\.

- log.p:

  logical; if `TRUE`, probabilities are returned on the log scale.

- n:

  number of random values to return (for `rbetabinom`).

## Value

`dbetabinom` gives the density, `pbetabinom` gives the distribution
function, and `rbetabinom` generates random samples.

## Details

This implementation of `dbetabinom` allows for automatic differentiation
with `RTMB`.

\$\$P(X = k;\\ n, a, b) = \binom{n}{k} \frac{B(k+a,\\ n-k+b)}{B(a,\\
b)}, \quad k = 0, 1, \ldots, n.\$\$

The distribution function has no closed form and is computed by summing
the probability mass function over \\0, \ldots, q\\. It is AD-compatible
in the parameters, while `q` and `size` must be numeric data. This is
also what one-step-ahead (OSA) residuals via
`RTMB::`[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
need, so these are supported, e.g. with `method = "cdf"` and
`discrete = TRUE`.

## Examples

``` r
set.seed(123)
x <- rbetabinom(1, 10, 2, 5)
d <- dbetabinom(x, 10, 2, 5)
p <- pbetabinom(x, 10, 2, 5)
```
