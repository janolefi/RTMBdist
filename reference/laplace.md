# Laplace distribution

Density, distribution function, quantile function, and random generation
for the Laplace distribution.

## Usage

``` r
dlaplace(x, mu = 0, b = 1, eps = NULL, log = FALSE)

plaplace(q, mu = 0, b = 1, lower.tail = TRUE, log.p = FALSE)

qlaplace(p, mu = 0, b = 1, lower.tail = TRUE, log.p = FALSE)

rlaplace(n, mu = 0, b = 1)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter

- b:

  scale parameter, must be positive.

- eps:

  optional smoothing parameter for `dlaplace` to smooth the absolute
  value function. See
  [`abs_smooth`](https://janolefi.github.io/RTMBdist/reference/abs_smooth.md)
  for details. It is recommended to set this to a small constant like
  `1e-6` for numerical optimisation.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dlaplace` gives the density, `plaplace` gives the distribution
function, `qlaplace` gives the quantile function, and `rlaplace`
generates random deviates.

## Details

This implementation of `dlaplace` allows for automatic differentiation
with `RTMB`.

## Examples

``` r
x <- rlaplace(1, 1, 1)
d <- dlaplace(x, 1, 1)
p <- plaplace(x, 1, 1)
q <- qlaplace(p, 1, 1)
```
