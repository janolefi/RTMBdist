# Sichel distribution

Probability mass function, distribution function and random generation
for the Sichel distribution.

## Usage

``` r
dsichel(x, mu, sigma, nu, log = FALSE)

psichel(q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE)

rsichel(n, mu, sigma, nu)
```

## Arguments

- x:

  vector of non-negative counts.

- mu:

  mean parameter, must be positive.

- sigma:

  dispersion parameter, must be positive.

- nu:

  shape parameter, any real number.

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

  number of random values to return (for `rsichel`).

## Value

`dsichel` gives the probability mass function, `psichel` gives the
distribution function, and `rsichel` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

The Sichel distribution is a mixed Poisson distribution whose mean
follows a generalised inverse Gaussian distribution. In the
parameterisation of the `SICHEL` family of the `gamlss.dist` package,
\\\mu\\ is the mean and the probability mass function is \$\$P(X = x) =
\frac{(\mu/c)^x\\ K\_{x + \nu}(\alpha)}{x!\\(\alpha\sigma)^{x + \nu}
K\_\nu(1/\sigma)}, \qquad \alpha^2 = \frac{1}{\sigma^2} +
\frac{2\mu}{c\\\sigma}, \quad c = \frac{K\_{\nu +
1}(1/\sigma)}{K\_\nu(1/\sigma)},\$\$ for \\x = 0, 1, 2, \ldots\\, where
\\K\_\nu\\ is the modified Bessel function of the second kind. The
variance is \$\$\mathrm{Var}(X) = \mu + \mu^2\left(\frac{2\sigma(\nu +
1)}{c} + \frac{1}{c^2} - 1\right).\$\$ For \\\nu = -1/2\\ this is the
[Poisson-inverse
Gaussian](https://janolefi.github.io/RTMBdist/reference/pig.md)
distribution, and as \\\sigma \to \infty\\ with \\\nu \> 0\\ it
approaches the [negative
binomial](https://janolefi.github.io/RTMBdist/reference/nbinom2.md).

The Bessel functions are evaluated with `RTMB`'s `besselK` for orders up
to 50. For larger orders, `besselK` overflows when its argument is
small, so its uniform asymptotic expansion for large orders is used
instead, which is accurate to about \\10^{-10}\\ on the log scale there.
Both are differentiable in the order, so \\\nu\\ can be estimated.

The distribution function has no closed form and is computed by summing
the probability mass function over \\0, \ldots, q\\, so its cost grows
with the largest `q`. It is AD-compatible in the parameters, while `q`
must be numeric data. This is also what one-step-ahead (OSA) residuals
via
`RTMB::`[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
need, so these are supported, e.g. with `method = "cdf"` and
`discrete = TRUE`.

## References

Sichel, H. S. (1982) Repeat-buying and the generalized inverse
Gaussian-Poisson distribution. Journal of the Royal Statistical Society,
Series C, 31, 193-204.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

Olver, F. W. J. et al. (eds.) NIST Digital Library of Mathematical
Functions, equation 10.41.4, https://dlmf.nist.gov/10.41.

## See also

[pig](https://janolefi.github.io/RTMBdist/reference/pig.md),
[nbinom2](https://janolefi.github.io/RTMBdist/reference/nbinom2.md),
[bnbinom2](https://janolefi.github.io/RTMBdist/reference/bnbinom2.md)

## Examples

``` r
set.seed(123)
x <- rsichel(5, mu = 4, sigma = 0.5, nu = 1)
d <- dsichel(x, mu = 4, sigma = 0.5, nu = 1)
p <- psichel(x, mu = 4, sigma = 0.5, nu = 1)
```
