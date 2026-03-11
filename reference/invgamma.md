# Inverse Gamma distribution

Density, distribution function, and random generation for the inverse
Gamma distribution.

## Usage

``` r
dinvgamma(x, shape, rate, scale = 1/rate, log = FALSE)

pinvgamma(q, shape, rate, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)

qinvgamma(p, shape, rate, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)

rinvgamma(n, shape, rate, scale = 1/rate)
```

## Arguments

- x, q:

  vector of quantiles, must be positive.

- shape, rate, scale:

  positive parameters of corresponding gamma distribution

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

`dinvgamma` gives the density, `pinvgamma` gives the distribution
function, `qinvgamma` gives the quantile function, and `rinvgamma`
generates random deviates.

## Details

This implementation of `dinvgamma`, `pinvgamma`, and `qinvgamma` allows
for automatic differentiation with `RTMB`.

If \\X \sim \Gamma(\alpha, \beta)\\, then \\1/X \sim
\text{InvGamma}(\alpha, \beta)\\.

## Examples

``` r
x <- rinvgamma(1, 1, 0.5)
d <- dinvgamma(x, 1, 0.5)
p <- pinvgamma(x, 1, 0.5)
q <- qinvgamma(p, 1, 0.5)
```
