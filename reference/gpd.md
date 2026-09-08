# Generalised Pareto distribution

Density, distribution function, quantile function, and random generation
for the generalised Pareto distribution (GPD).

## Usage

``` r
dgpd(x, mu = 0, sigma = 1, xi = 0, log = FALSE)

pgpd(q, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE)

qgpd(p, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE)

rgpd(n, mu = 0, sigma = 1, xi = 0)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter, the threshold below which the density is zero.

- sigma:

  scale parameter, must be positive.

- xi:

  shape parameter (real).

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dgpd` gives the density, `pgpd` gives the distribution function, `qgpd`
gives the quantile function, and `rgpd` generates random deviates.

## Details

`dgpd` and `pgpd` allow for automatic differentiation with `RTMB`.

With \\z = (x - \mu) / \sigma\\ the survival function is \$\$1 -
F(x;\\\mu,\sigma,\xi) = \begin{cases} (1 + \xi z)^{-1/\xi} & \xi \neq 0,
\\ e^{-z} & \xi = 0, \end{cases}\$\$ for \\x \ge \mu\\, and the density
is \$\$f(x;\\\mu,\sigma,\xi) = \frac{1}{\sigma} (1 + \xi z)^{-1/\xi -
1}.\$\$

This is the limiting distribution of exceedances over a high threshold,
so \\\mu\\ is usually a fixed threshold rather than an estimated
parameter. The support is \\x \ge \mu\\ for \\\xi \ge 0\\ and \\\mu \le
x \le \mu - \sigma/\xi\\ for \\\xi \< 0\\. At \\\xi = 0\\ the
distribution is exponential with rate \\1/\sigma\\, and at \\\xi \> 0\\
with \\\mu = \sigma/\xi\\ it is the
[Pareto](https://janolefi.github.io/RTMBdist/reference/pareto.md)
distribution with \\\mu = 1/\xi\\.

The three cases are covered by one expression, so no branch on the sign
or the value of \\\xi\\ is needed. In particular the derivative with
respect to \\\xi\\ is exact at \\\xi = 0\\, which is the usual starting
value when the shape is estimated.

The threshold itself belongs to the support: `dgpd(mu, mu, sigma, xi)`
is \\1/\sigma\\, as `stats::dexp(0, rate)` is `rate`. The `VGAM`, `evd`
and `extraDistr` implementations return zero there instead.

## References

Coles, S. (2001) An Introduction to Statistical Modeling of Extreme
Values, Springer, doi:10.1007/978-1-4471-3675-0.

Pickands, J. (1975) Statistical inference using extreme order
statistics. The Annals of Statistics, 3, 119-131.

## See also

[gev](https://janolefi.github.io/RTMBdist/reference/gev.md),
[pareto](https://janolefi.github.io/RTMBdist/reference/pareto.md),
[frechet](https://janolefi.github.io/RTMBdist/reference/frechet.md)

## Examples

``` r
set.seed(123)
x <- rgpd(5, mu = 0, sigma = 1, xi = 0.3)
d <- dgpd(x, mu = 0, sigma = 1, xi = 0.3)
p <- pgpd(x, mu = 0, sigma = 1, xi = 0.3)
q <- qgpd(p, mu = 0, sigma = 1, xi = 0.3)
```
