# Skew normal distribution

Density, distribution function, quantile function, and random generation
for the skew normal distribution.

## Usage

``` r
dskewnorm(x, xi = 0, omega = 1, alpha = 0, log = FALSE)

pskewnorm(q, xi = 0, omega = 1, alpha = 0, ...)

qskewnorm(p, xi = 0, omega = 1, alpha = 0, ...)

rskewnorm(n, xi = 0, omega = 1, alpha = 0)
```

## Arguments

- x, q:

  vector of quantiles

- xi:

  location parameter

- omega:

  scale parameter, must be positive.

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

`dskewnorm` gives the density, `pskewnorm` gives the distribution
function, `qskewnorm` gives the quantile function, and `rskewnorm`
generates random deviates.

## Details

This implementation of `dskewnorm` allows for automatic differentiation
with `RTMB` while the other functions are imported from the `sn`
package. See `sn::`[`dsn`](https://rdrr.io/pkg/sn/man/dsn.html) for more
details.

## See also

[skewnorm2](https://janolefi.github.io/RTMBdist/reference/skewnorm2.md),
[skewt](https://janolefi.github.io/RTMBdist/reference/skewt.md),
[skewt2](https://janolefi.github.io/RTMBdist/reference/skewt2.md)

## Examples

``` r
# alpha is skew parameter
x <- rskewnorm(1, alpha = 1)
d <- dskewnorm(x, alpha = 1)
p <- pskewnorm(x, alpha = 1)
q <- qskewnorm(p, alpha = 1)
```
