# Wrapped normal distribution

Density, distribution function, and random generation for the wrapped
normal distribution.

## Usage

``` r
dwrpnorm(x, mu = 0, sigma, log = FALSE)

pwrpnorm(q, mu = 0, sigma, from = NULL, lower.tail = TRUE, log.p = FALSE)

rwrpnorm(n, mu = 0, sigma)
```

## Arguments

- x, q:

  vector of angles in radians.

- mu:

  mean direction in radians.

- sigma:

  standard deviation of the normal distribution before wrapping, must be
  positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- from:

  origin, in radians, at which the circle is cut open for the
  distribution function. If `NULL` (default), it is set to `mu - pi`.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dwrpnorm` gives the density, `pwrpnorm` gives the distribution
function, and `rwrpnorm` generates random deviates in \\\[-\pi, \pi)\\.

## Details

This implementation allows for automatic differentiation with `RTMB`.

The wrapped normal is the distribution of \\Y \bmod 2\pi\\ for a normal
\\Y\\ with mean \\\mu\\ and standard deviation \\\sigma\\, with density
\$\$f(x;\\\mu,\sigma) = \sum\_{k=-\infty}^{\infty}
\frac{1}{\sigma}\\\phi\\\left(\frac{x - \mu + 2\pi k}{\sigma}\right) =
\frac{1}{2\pi}\left(1 + 2\sum\_{k=1}^{\infty} e^{-k^2\sigma^2/2}
\cos(k(x - \mu))\right).\$\$ For \\\sigma \< 2\\ the first sum is used
with \\\|k\| \le 3\\, otherwise the second with \\k \le 4\\; the terms
left out are below \\10^{-20}\\.

The distribution function cuts the circle open at \\\mu - \pi\\ by
default, or at `from`, exactly as
[`pwrpcauchy`](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md)
does, where the choice of origin and the one-step-ahead (OSA) residuals
are explained. As there, OSA residuals via
`RTMB::`[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
with `method = "cdf"` cut the circle at the fixed origin \\-\pi\\.

## See also

[wrpcauchy](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md),
[vm](https://janolefi.github.io/RTMBdist/reference/vm.md),
[card](https://janolefi.github.io/RTMBdist/reference/card.md);
[`cjw()`](https://janolefi.github.io/RTMBdist/reference/cjw.md) for
circular-linear copulas.

## Examples

``` r
set.seed(1)
x <- rwrpnorm(10, 0, 0.8)
d <- dwrpnorm(x, 0, 0.8)
p <- pwrpnorm(x, 0, 0.8)
```
