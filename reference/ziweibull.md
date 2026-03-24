# Zero-inflated Weibull distribution

Density, distribution function, and random generation for the
zero-inflated Weibull distribution.

## Usage

``` r
dziweibull(x, shape, scale, zeroprob = 0, log = FALSE)

pziweibull(q, shape, scale, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rziweibull(n, shape, scale, zeroprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- shape:

  positive shape parameter

- scale:

  positive scale parameter

- zeroprob:

  zero-inflation probability between 0 and 1.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return

## Value

`dziweibull` gives the density, `pziweibull` gives the distribution
function, and `rziweibull` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

## Examples

``` r
x <- rziweibull(1, 1, 1, 0.5)
d <- dziweibull(x, 1, 1, 0.5)
p <- pziweibull(x, 1, 1, 0.5)
```
