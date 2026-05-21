# Exponentially modified Gaussian distribution

Density, distribution function, quantile function, and random generation
for the exponentially modified Gaussian distribution.

## Usage

``` r
dexgauss(x, mu = 0, sigma = 1, lambda = 1, log = FALSE)

pexgauss(q, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE)

qexgauss(p, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE)

rexgauss(n, mu = 0, sigma = 1, lambda = 1)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  mean parameter of the Gaussian part

- sigma:

  standard deviation parameter of the Gaussian part, must be positive.

- lambda:

  rate parameter of the exponential part, must be positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dexgauss` gives the density, `pexgauss` gives the distribution
function, `qexgauss` gives the quantile function, and `rexgauss`
generates random deviates.

## Details

This implementation of `dexgauss` and `pexgauss` allows for automatic
differentiation with `RTMB`. `qexgauss` and `rexgauss` are
reparameterised imports from
`gamlss.dist::`[`exGAUS`](https://rdrr.io/pkg/gamlss.dist/man/exGAUS.html).

If \\X \sim N(\mu, \sigma^2)\\ and \\Y \sim \text{Exp}(\lambda)\\, then
\\Z = X + Y\\ follows the exponentially modified Gaussian distribution
with parameters \\\mu\\, \\\sigma\\, and \\\lambda\\.

The density is \$\$f(x;\\\mu,\sigma,\lambda) = \lambda
\exp\\\Bigl(\lambda\mu + \tfrac{\lambda^2\sigma^2}{2} - \lambda
x\Bigr)\\ \Phi\\\left(\frac{x - \mu -
\lambda\sigma^2}{\sigma}\right),\$\$ where \\\Phi\\ is the standard
normal CDF.

## References

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## Examples

``` r
x <- rexgauss(1, 1, 2, 2)
d <- dexgauss(x, 1, 2, 2)
p <- pexgauss(x, 1, 2, 2)
q <- qexgauss(p, 1, 2, 2)
```
