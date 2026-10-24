# Poisson-inverse Gaussian distribution

Probability mass function, distribution function and random generation
for the Poisson-inverse Gaussian distribution.

## Usage

``` r
dpig(x, mu, sigma, log = FALSE)

ppig(q, mu, sigma, lower.tail = TRUE, log.p = FALSE)

rpig(n, mu, sigma)
```

## Arguments

- x:

  vector of non-negative counts.

- mu:

  mean parameter, must be positive.

- sigma:

  dispersion parameter, must be positive.

- log:

  logical; if `TRUE`, probabilities are returned on the log scale.

- q:

  vector of quantiles.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le q\]\\,
  otherwise \\P\[X \> q\]\\.

- log.p:

  logical; if `TRUE`, probabilities are returned on the log scale.

- n:

  number of random values to return (for `rpig`).

## Value

`dpig` gives the probability mass function, `ppig` gives the
distribution function, and `rpig` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

The Poisson-inverse Gaussian is a mixed Poisson distribution whose mean
is drawn from an inverse Gaussian distribution with mean \\\mu\\ and
variance \\\sigma \mu^2\\. Its probability mass function is \$\$P(X = x)
= \left(\frac{2\alpha}{\pi}\right)^{1/2} \frac{\mu^x e^{1/\sigma}
K\_{x - 1/2}(\alpha)}{(\alpha\sigma)^x\\ x!}, \qquad \alpha^2 =
\frac{1}{\sigma^2} + \frac{2\mu}{\sigma}, \quad x = 0, 1, 2, \ldots\$\$
where \\K\_\nu\\ is the modified Bessel function of the second kind. The
parameterisation follows the `PIG` family of the `gamlss.dist` package,
with \$\$E(X) = \mu, \qquad \mathrm{Var}(X) = \mu + \sigma \mu^2.\$\$
This is the same mean-variance relationship as the [negative
binomial](https://janolefi.github.io/RTMBdist/reference/nbinom2.md) with
`size = 1 / sigma`, but the Poisson-inverse Gaussian has a heavier right
tail and puts more mass at zero for the same mean and variance. It
approaches the Poisson as \\\sigma \to 0\\.

It is the
[Sichel](https://janolefi.github.io/RTMBdist/reference/sichel.md)
distribution with \\\nu = -1/2\\, and is computed as such; see there for
how the Bessel functions are evaluated. The distribution function is
summed from the probability mass function, so one-step-ahead (OSA)
residuals via
`RTMB::`[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
are supported, e.g. with `method = "cdf"` and `discrete = TRUE`.

## References

Dean, C., Lawless, J. F., and Willmot, G. E. (1989) A mixed
Poisson-inverse-Gaussian regression model. Canadian Journal of
Statistics, 17, 171-181.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[sichel](https://janolefi.github.io/RTMBdist/reference/sichel.md),
[nbinom2](https://janolefi.github.io/RTMBdist/reference/nbinom2.md),
[invgauss](https://janolefi.github.io/RTMBdist/reference/invgauss.md)

## Examples

``` r
set.seed(123)
x <- rpig(5, mu = 4, sigma = 0.5)
d <- dpig(x, mu = 4, sigma = 0.5)
p <- ppig(x, mu = 4, sigma = 0.5)
```
