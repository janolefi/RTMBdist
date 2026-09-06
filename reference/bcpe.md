# Box-Cox Power Exponential distribution (BCPE)

Density, distribution function, quantile function, and random generation
for the Box-Cox Power Exponential distribution.

## Usage

``` r
dbcpe(x, mu = 5, sigma = 0.1, nu = 1, tau = 2, log = FALSE)

pbcpe(q, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)

qbcpe(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)

rbcpe(n, mu = 5, sigma = 0.1, nu = 1, tau = 2)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter, must be positive.

- sigma:

  scale parameter, must be positive.

- nu:

  vector of `nu` parameter values.

- tau:

  vector of `tau` parameter values, must be positive.

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

`dbcpe` gives the density, `pbcpe` gives the distribution function,
`qbcpe` gives the quantile function, and `rbcpe` generates random
deviates.

## Details

`dbcpe` and `pbcpe` allow for automatic differentiation with `RTMB`. The
parameterisation follows the `BCPE` family of the `gamlss.dist` package.

The density is \$\$f(x; \mu, \sigma, \nu, \tau) =
\frac{x^{\nu-1}}{\mu^{\nu} \sigma}
\frac{f_T(z;\tau)}{F_T\\\left(1/(\sigma\|\nu\|);\tau\right)}, \quad x \>
0,\$\$ where \\z = \[(x/\mu)^\nu - 1\]/(\nu\sigma)\\ for \\\nu \neq 0\\
and \\z = \log(x/\mu)/\sigma\\ for \\\nu = 0\\, and \\f_T(\cdot;\tau)\\
and \\F_T(\cdot;\tau)\\ are the PDF and CDF of the power exponential
(PE) distribution with shape \\\tau\\.

## References

Rigby, R. A. and Stasinopoulos, D. M. (2004) Smooth centile curves for
skew and kurtotic data modelled using the Box-Cox Power Exponential
distribution. Statistics in Medicine, 23, 3053-3076.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[bccg](https://janolefi.github.io/RTMBdist/reference/bccg.md),
[bct](https://janolefi.github.io/RTMBdist/reference/bct.md),
[powerexp](https://janolefi.github.io/RTMBdist/reference/powerexp.md)

## Examples

``` r
x <- rbcpe(1, mu = 5, sigma = 0.1, nu = 1, tau = 1)
d <- dbcpe(x, mu = 5, sigma = 0.1, nu = 1, tau = 1)
p <- pbcpe(x, mu = 5, sigma = 0.1, nu = 1, tau = 1)
q <- qbcpe(p, mu = 5, sigma = 0.1, nu = 1, tau = 1)
```
