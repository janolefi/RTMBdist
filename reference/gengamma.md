# Generalised Gamma distribution (GG)

Density, distribution function, quantile function, and random generation
for the generalised Gamma distribution.

## Usage

``` r
dgengamma(x, mu = 1, sigma = 0.5, nu = 1, log = FALSE)

pgengamma(q, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)

qgengamma(p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, log.p = FALSE)

rgengamma(n, mu = 1, sigma = 0.5, nu = 1)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter, must be positive.

- sigma:

  scale parameter, must be positive.

- nu:

  skewness parameter (real).

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

`dgengamma` gives the density, `pgengamma` gives the distribution
function, `qgengamma` gives the quantile function, and `rgengamma`
generates random deviates.

## Details

This implementation of `dgengamma`, `pgengamma`, and `qgengamma` allows
for automatic differentiation with `RTMB`.

## References

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## Examples

``` r
x <- rgengamma(5, mu = 4, sigma = 0.5, nu = 0.5)
d <- dgengamma(x, mu = 4, sigma = 0.5, nu = 0.5)
p <- pgengamma(x, mu = 4, sigma = 0.5, nu = 0.5)
q <- qgengamma(p, mu = 4, sigma = 0.5, nu = 0.5)
```
