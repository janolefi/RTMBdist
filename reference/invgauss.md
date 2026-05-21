# Inverse Gaussian distribution

Density, distribution function, and random generation for the inverse
Gaussian distribution.

## Usage

``` r
dinvgauss(x, mean = 1, shape = 1, log = FALSE)

pinvgauss(q, mean = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)

qinvgauss(p, mean = 1, shape = 1, lower.tail = TRUE, log.p = FALSE, ...)

rinvgauss(n, mean = 1, shape = 1)
```

## Arguments

- x, q:

  vector of quantiles, must be positive.

- mean:

  location parameter

- shape:

  shape parameter, must be positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- ...:

  additional parameter passed to
  [`statmod::qinvgauss`](https://rdrr.io/pkg/statmod/man/invgauss.html)
  for numerical evaluation of the quantile function.

- n:

  number of random values to return

## Value

`dinvgauss` gives the density, `pinvgauss` gives the distribution
function, `qinvgauss` gives the quantile function, and `rinvgauss`
generates random deviates.

## Details

This implementation of `dinvgauss` allows for automatic differentiation
with `RTMB`. `qinvgauss` and `rinvgauss` are imported from the `statmod`
package.

\$\$f(x;\\\mu,\lambda) = \sqrt{\frac{\lambda}{2\pi x^3}}
\exp\\\left(-\frac{\lambda(x-\mu)^2}{2\mu^2 x}\right), \quad x \> 0.\$\$

## Examples

``` r
x <- rinvgauss(1, 1, 0.5)
d <- dinvgauss(x, 1, 0.5)
p <- pinvgauss(x, 1, 0.5)
q <- qinvgauss(p, 1, 0.5)
```
