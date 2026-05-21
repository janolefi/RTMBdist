# Reparameterised one-inflated beta distribution

Density, distribution function, and random generation for the
one-inflated beta distribution reparameterised in terms of mean and
concentration.

## Usage

``` r
doibeta2(x, mu, phi, oneprob = 0, log = FALSE)

poibeta2(q, mu, phi, oneprob = 0, lower.tail = TRUE, log.p = FALSE)

roibeta2(n, mu, phi, oneprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  mean parameter, must be in the interval from 0 to 1.

- phi:

  concentration parameter, must be positive.

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

`doibeta2` gives the density, `poibeta2` gives the distribution
function, and `roibeta2` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

## Examples

``` r
set.seed(123)
x <- roibeta2(1, 0.6, 2, 0.5)
d <- doibeta2(x, 0.6, 2, 0.5)
p <- poibeta2(x, 0.6, 2, 0.5)
```
