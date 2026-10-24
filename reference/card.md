# Cardioid distribution

Density, distribution function, and random generation for the cardioid
distribution.

## Usage

``` r
dcard(x, mu = 0, rho, log = FALSE)

pcard(q, mu = 0, rho, from = NULL, lower.tail = TRUE, log.p = FALSE)

rcard(n, mu = 0, rho)
```

## Arguments

- x, q:

  vector of angles in radians.

- mu:

  mean direction in radians.

- rho:

  concentration parameter, must be in the interval from 0 to 1/2.

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

`dcard` gives the density, `pcard` gives the distribution function, and
`rcard` generates random deviates in \\\[-\pi, \pi)\\.

## Details

This implementation allows for automatic differentiation with `RTMB`.

The cardioid distribution has density \$\$f(x;\\\mu,\rho) = \frac{1 +
2\rho\cos(x - \mu)}{2\pi}, \qquad 0 \le \rho \le \frac{1}{2},\$\$ and is
the circular uniform for \\\rho = 0\\. With the circle cut open at
\\\mu - \pi\\, its distribution function is \$\$F(q) = \frac{q - \mu +
\pi}{2\pi} + \frac{\rho}{\pi}\sin(q - \mu), \quad \mu - \pi \le q \<
\mu + \pi.\$\$ The origin can be changed with `from`, exactly as for
[`pwrpcauchy`](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md),
where the choice of origin and the one-step-ahead (OSA) residuals are
explained. As there, OSA residuals via
`RTMB::`[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
with `method = "cdf"` cut the circle at the fixed origin \\-\pi\\.

## References

Jammalamadaka, S. R. and SenGupta, A. (2001) Topics in Circular
Statistics. World Scientific.

## See also

[wrpcauchy](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md),
[wrpnorm](https://janolefi.github.io/RTMBdist/reference/wrpnorm.md),
[vm](https://janolefi.github.io/RTMBdist/reference/vm.md)

## Examples

``` r
set.seed(1)
x <- rcard(10, 0, 0.3)
d <- dcard(x, 0, 0.3)
p <- pcard(x, 0, 0.3)
```
