# One-inflated beta distribution

Density, distribution function, and random generation for the
one-inflated beta distribution.

## Usage

``` r
doibeta(x, shape1, shape2, oneprob = 0, log = FALSE)

poibeta(q, shape1, shape2, oneprob = 0, lower.tail = TRUE, log.p = FALSE)

roibeta(n, shape1, shape2, oneprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- shape1, shape2:

  non-negative shape parameters of the beta distribution

- oneprob:

  one-inflation probability between 0 and 1.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`doibeta` gives the density, `poibeta` gives the distribution function,
and `roibeta` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

## Examples

``` r
set.seed(123)
x <- roibeta(1, 2, 2, 0.5)
d <- doibeta(x, 2, 2, 0.5)
p <- poibeta(x, 2, 2, 0.5)
```
