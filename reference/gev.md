# Generalised extreme value distribution

Density, distribution function, quantile function, and random generation
for the generalised extreme value (GEV) distribution.

## Usage

``` r
dgev(x, mu = 0, sigma = 1, xi = 0, log = FALSE)

pgev(q, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE)

qgev(p, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE)

rgev(n, mu = 0, sigma = 1, xi = 0)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter

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

`dgev` gives the density, `pgev` gives the distribution function, `qgev`
gives the quantile function, and `rgev` generates random deviates.

## Details

`dgev` and `pgev` allow for automatic differentiation with `RTMB`.

With \\z = (x - \mu) / \sigma\\ the distribution function is
\$\$F(x;\\\mu,\sigma,\xi) = \exp\bigl(-t(x)\bigr), \qquad t(x) =
\begin{cases} (1 + \xi z)^{-1/\xi} & \xi \neq 0, \\ e^{-z} & \xi = 0,
\end{cases}\$\$ and the density is \\f(x) = t(x)^{\xi + 1} e^{-t(x)} /
\sigma\\.

The shape \\\xi\\ determines the tail and with it the support: for \\\xi
\> 0\\ (Frechet case) the distribution is heavy-tailed on \\x \> \mu -
\sigma/\xi\\, for \\\xi \< 0\\ (Weibull case) it is bounded above by
\\\mu - \sigma/\xi\\, and \\\xi = 0\\ is the
[Gumbel](https://janolefi.github.io/RTMBdist/reference/gumbel.md)
distribution on the whole real line.

The three cases are covered by one expression, so no branch on the sign
or the value of \\\xi\\ is needed. In particular the derivative with
respect to \\\xi\\ is exact at \\\xi = 0\\, which is the usual starting
value when the shape is estimated.

## References

Coles, S. (2001) An Introduction to Statistical Modeling of Extreme
Values, Springer, doi:10.1007/978-1-4471-3675-0.

Jenkinson, A. F. (1955) The frequency distribution of the annual maximum
(or minimum) values of meteorological elements. Quarterly Journal of the
Royal Meteorological Society, 81, 158-171.

## See also

[gumbel](https://janolefi.github.io/RTMBdist/reference/gumbel.md),
[gpd](https://janolefi.github.io/RTMBdist/reference/gpd.md),
[frechet](https://janolefi.github.io/RTMBdist/reference/frechet.md)

## Examples

``` r
set.seed(123)
x <- rgev(5, mu = 0, sigma = 1, xi = 0.2)
d <- dgev(x, mu = 0, sigma = 1, xi = 0.2)
p <- pgev(x, mu = 0, sigma = 1, xi = 0.2)
q <- qgev(p, mu = 0, sigma = 1, xi = 0.2)
```
