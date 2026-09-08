# Yule-Simon distribution

Probability mass function, distribution function, and random generation
for the Yule-Simon distribution.

## Usage

``` r
dyules(x, shape = 1, log = FALSE)

pyules(q, shape = 1, lower.tail = TRUE, log.p = FALSE)

ryules(n, shape = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- shape:

  positive shape parameter \\\rho\\.

- log, log.p:

  logical; if `TRUE`, probabilities are returned as \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- n:

  number of random values to return (for `ryules`).

## Value

`dyules` gives the probability mass function, `pyules` gives the
distribution function, and `ryules` generates random deviates.

## Details

`dyules` and `pyules` allow for automatic differentiation with `RTMB`.

\$\$P(X = k;\\ \rho) = \rho\\ B(k,\\ \rho + 1), \quad k = 1, 2,
\ldots\$\$

The Yule-Simon distribution is the classical long-tailed frequency law,
used for word counts, city sizes, citation counts and species-per-genus
data. Its tail is a power law, \\P(X = k) \sim \rho\\\Gamma(\rho + 1)
k^{-(\rho + 1)}\\, so the mean \\\rho/(\rho - 1)\\ exists only for
\\\rho \> 1\\ and the variance \\\rho^2 / \\(\rho - 1)^2 (\rho - 2)\\\\
only for \\\rho \> 2\\.

It is the [beta-negative
binomial](https://janolefi.github.io/RTMBdist/reference/bnbinom.md) with
`size = 1` and `shape2 = 1`, shifted to start at one, and equivalently
the [Waring](https://janolefi.github.io/RTMBdist/reference/waring.md)
distribution with `sigma = mu`, shifted the same way. That special case
is what gives it a closed-form distribution function, which the
beta-negative binomial does not have in general.

The support here starts at one, as in `VGAM`. The `YULE` family of
`gamlss.dist` instead starts at zero and is parameterised by its mean
\\\mu\\, which corresponds to \\\rho = (\mu + 1)/\mu\\; use
`dwaring(x, mu, mu)` for that version.

## References

Simon, H. A. (1955) On a class of skew distribution functions.
Biometrika, 42, 425-440, doi:10.1093/biomet/42.3-4.425.

## See also

[waring](https://janolefi.github.io/RTMBdist/reference/waring.md),
[bnbinom](https://janolefi.github.io/RTMBdist/reference/bnbinom.md),
[bnbinom2](https://janolefi.github.io/RTMBdist/reference/bnbinom2.md)

## Examples

``` r
set.seed(123)
x <- ryules(5, shape = 2)
d <- dyules(x, shape = 2)
p <- pyules(x, shape = 2)
```
