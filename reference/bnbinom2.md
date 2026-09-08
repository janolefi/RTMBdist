# Reparameterised beta-negative binomial distribution

Probability mass function and random generation for the beta-negative
binomial distribution reparameterised in terms of its mean.

## Usage

``` r
dbnbinom2(x, mu, sigma, nu, log = FALSE)

rbnbinom2(n, mu, sigma, nu)
```

## Arguments

- x:

  vector of non-negative counts.

- mu:

  mean parameter, must be positive.

- sigma:

  dispersion parameter, must be positive. The variance is finite only
  for `sigma < 1`.

- nu:

  dispersion parameter, must be positive.

- log:

  logical; if `TRUE`, probabilities are returned on the log scale.

- n:

  number of random values to return (for `rbnbinom2`).

## Value

`dbnbinom2` gives the probability mass function and `rbnbinom2`
generates random deviates.

## Details

`dbnbinom2` allows for automatic differentiation with `RTMB`. The
parameterisation follows the `BNB` family of the `gamlss.dist` package,
in which \\\mu\\ is exactly the mean.

Writing \\r\\, \\a\\ and \\b\\ for the arguments `size`, `shape1` and
`shape2` of
[`dbnbinom`](https://janolefi.github.io/RTMBdist/reference/bnbinom.md),
the reparameterisation is \$\$a = \frac{1}{\sigma} + 1, \qquad b =
\frac{\mu\nu}{\sigma}, \qquad r = \frac{1}{\nu}.\$\$

This gives \\E(X) = \mu\\ for every admissible \\\sigma\\ and \\\nu\\,
where the original parameterisation needs \\a \> 1\\ for a mean to exist
at all. The variance is \$\$\mathrm{Var}(X) = \frac{\mu (\sigma +
\nu)(\mu\nu + 1)}{\nu (1 - \sigma)},\$\$ which is finite for \\\sigma \<
1\\ and increases without bound as \\\sigma \to 1\\. Both \\\sigma\\ and
\\\nu\\ add dispersion beyond the [negative
binomial](https://janolefi.github.io/RTMBdist/reference/nbinom2.md),
which is recovered as \\\sigma \to 0\\.

Because the mean is pinned to a single parameter, this parameterisation
is the more stable of the two to estimate in; see
[`bnbinom`](https://janolefi.github.io/RTMBdist/reference/bnbinom.md)
for the identifiability problem it avoids.

There is no distribution function, since the beta-negative binomial
distribution function has no closed form and the support is unbounded.
One-step-ahead residuals are therefore not available.

## References

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[bnbinom](https://janolefi.github.io/RTMBdist/reference/bnbinom.md),
[nbinom2](https://janolefi.github.io/RTMBdist/reference/nbinom2.md),
[betabinom](https://janolefi.github.io/RTMBdist/reference/betabinom.md)

## Examples

``` r
set.seed(123)
x <- rbnbinom2(5, mu = 4, sigma = 0.4, nu = 0.5)
d <- dbnbinom2(x, mu = 4, sigma = 0.4, nu = 0.5)
```
