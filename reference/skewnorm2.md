# Reparameterised skew normal distribution

Density, distribution function, quantile function and random generation
for the skew normal distribution reparameterised in terms of mean,
standard deviation and skew magnitude

## Usage

``` r
dskewnorm2(x, mean = 0, sd = 1, alpha = 0, log = FALSE)

pskewnorm2(q, mean = 0, sd = 1, alpha = 0, ...)

qskewnorm2(p, mean = 0, sd = 1, alpha = 0, ...)

rskewnorm2(n, mean = 0, sd = 1, alpha = 0)
```

## Arguments

- x, q:

  vector of quantiles

- mean:

  mean parameter

- sd:

  standard deviation, must be positive.

- alpha:

  skewness parameter, +/- `Inf` is allowed.

- log:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- ...:

  additional parameters to be passed to the `sn` package functions for
  `pskewnorm` and `qskewnorm`.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dskewnorm2` gives the density, `pskewnorm2` gives the distribution
function, `qskewnorm2` gives the quantile function, and `rskewnorm2`
generates random deviates.

## Details

This implementation of `dskewnorm2` allows for automatic differentiation
with `RTMB` while the other functions are imported from the `sn`
package.

## See also

[skewnorm](https://janolefi.github.io/RTMBdist/reference/skewnorm.md),
[skewt](https://janolefi.github.io/RTMBdist/reference/skewt.md),
[skewt2](https://janolefi.github.io/RTMBdist/reference/skewt2.md)

## Examples

``` r
# alpha is skew parameter
x <- rskewnorm2(1, alpha = 1)
d <- dskewnorm2(x, alpha = 1)
p <- pskewnorm2(x, alpha = 1)
q <- qskewnorm2(p, alpha = 1)
```
