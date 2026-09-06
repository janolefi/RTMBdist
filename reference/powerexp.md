# Power Exponential distribution (PE and PE2)

Density, distribution function, quantile function, and random generation
for the Power Exponential distribution (two versions).

## Usage

``` r
dpowerexp(x, mu = 0, sigma = 1, nu = 2, log = FALSE)

ppowerexp(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)

qpowerexp(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)

rpowerexp(n, mu = 0, sigma = 1, nu = 2)

dpowerexp2(x, mu = 0, sigma = 1, nu = 2, log = FALSE)

ppowerexp2(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)

qpowerexp2(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)

rpowerexp2(n, mu = 0, sigma = 1, nu = 2)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter

- sigma:

  scale parameter, must be positive

- nu:

  shape parameter (real)

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dpowerexp` gives the density, `ppowerexp` gives the distribution
function, `qpowerexp` gives the quantile function, and `rpowerexp`
generates random deviates.

## Details

The densities and distribution functions allow for automatic
differentiation with `RTMB`.

For `powerexp`, `mu` is the mean and `sigma` is the standard deviation
while this does not hold for `powerexp2`.

The parameterisation follows the `PE` and `PE2` families of the
`gamlss.dist` package.

For `powerexp` (PE), \\\sigma\\ is the standard deviation; the density
is \$\$f(x;\\\mu,\sigma,\nu) = \frac{\nu}{2c\\\sigma\\\Gamma(1/\nu)}
\exp\\\left(-\tfrac{1}{2}\left\|\frac{x-\mu}{c\sigma}\right\|^\nu\right),\$\$
where \\c = \[2^{-2/\nu}\Gamma(1/\nu)/\Gamma(3/\nu)\]^{1/2}\\.

For `powerexp2` (PE2), \\\sigma\\ is a scale parameter; the density is
\$\$f(x;\\\mu,\sigma,\nu) = \frac{\nu}{2\sigma\\\Gamma(1/\nu)}
\exp\\\left(-\left\|\frac{x-\mu}{\sigma}\right\|^\nu\right).\$\$

## References

Nelson, D. B. (1991) Conditional heteroskedasticity in asset returns: a
new approach. Econometrica, 59, 347-370.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[bcpe](https://janolefi.github.io/RTMBdist/reference/bcpe.md),
[jsu](https://janolefi.github.io/RTMBdist/reference/jsu.md),
[laplace](https://janolefi.github.io/RTMBdist/reference/laplace.md)

## Examples

``` r
# PE
x <- rpowerexp(1, mu = 0, sigma = 1, nu = 2)
d <- dpowerexp(x, mu = 0, sigma = 1, nu = 2)
p <- ppowerexp(x, mu = 0, sigma = 1, nu = 2)
q <- qpowerexp(p, mu = 0, sigma = 1, nu = 2)

# PE2
x <- rpowerexp2(1, mu = 0, sigma = 1, nu = 2)
d <- dpowerexp2(x, mu = 0, sigma = 1, nu = 2)
p <- ppowerexp2(x, mu = 0, sigma = 1, nu = 2)
q <- qpowerexp2(p, mu = 0, sigma = 1, nu = 2)
```
