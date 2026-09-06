# Box–Cox Cole and Green distribution (BCCG)

Density, distribution function, quantile function, and random generation
for the Box–Cox Cole and Green distribution.

## Usage

``` r
dbccg(x, mu = 1, sigma = 0.1, nu = 1, log = FALSE)

pbccg(q, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE)

qbccg(p, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE)

rbccg(n, mu = 1, sigma = 0.1, nu = 1)
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

`dbccg` gives the density, `pbccg` gives the distribution function,
`qbccg` gives the quantile function, and `rbccg` generates random
deviates.

## Details

`dbccg` and `pbccg` allow for automatic differentiation with `RTMB`. The
parameterisation follows the `BCCG` family of the `gamlss.dist` package.

The density is \$\$f(x; \mu, \sigma, \nu) = \frac{x^{\nu-1}}{\mu^{\nu}
\sigma \sqrt{2\pi}} \exp\\\left(-\tfrac{z^2}{2}\right)
\Bigl\[\Phi\\\left(\tfrac{1}{\sigma\|\nu\|}\right)\Bigr\]^{-1}, \quad x
\> 0,\$\$ where \\z = \[(x/\mu)^\nu - 1\]/(\nu\sigma)\\ for \\\nu \neq
0\\ and \\z = \log(x/\mu)/\sigma\\ for \\\nu = 0\\, and \\\Phi\\ is the
standard normal CDF.

## References

Cole, T. J. and Green, P. J. (1992) Smoothing reference centile curves:
the LMS method and penalized likelihood. Statistics in Medicine, 11,
1305-1319.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[bct](https://janolefi.github.io/RTMBdist/reference/bct.md),
[bcpe](https://janolefi.github.io/RTMBdist/reference/bcpe.md),
[gengamma](https://janolefi.github.io/RTMBdist/reference/gengamma.md)

## Examples

``` r
x <- rbccg(5, mu = 10, sigma = 0.2, nu = 0.5)
d <- dbccg(x, mu = 10, sigma = 0.2, nu = 0.5)
p <- pbccg(x, mu = 10, sigma = 0.2, nu = 0.5)
q <- qbccg(p, mu = 10, sigma = 0.2, nu = 0.5)
```
