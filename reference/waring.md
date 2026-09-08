# Waring distribution

Probability mass function, distribution function, and random generation
for the Waring distribution.

## Usage

``` r
dwaring(x, mu = 2, sigma = 2, log = FALSE)

pwaring(q, mu = 2, sigma = 2, lower.tail = TRUE, log.p = FALSE)

rwaring(n, mu = 2, sigma = 2)
```

## Arguments

- x, q:

  vector of non-negative counts.

- mu:

  mean parameter, must be positive.

- sigma:

  dispersion parameter, must be positive. The variance is finite only
  for `sigma < 1`.

- log, log.p:

  logical; if `TRUE`, probabilities are returned as \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- n:

  number of random values to return (for `rwaring`).

## Value

`dwaring` gives the probability mass function, `pwaring` gives the
distribution function, and `rwaring` generates random deviates.

## Details

`dwaring` and `pwaring` allow for automatic differentiation with `RTMB`.
The parameterisation follows the `WARING` family of the `gamlss.dist`
package, in which \\\mu\\ is exactly the mean.

\$\$P(X = k;\\ \mu, \sigma) = \frac{B\bigl(k + \tfrac{\mu}{\sigma},\\
\tfrac{1}{\sigma} + 2\bigr)}{B\bigl(\tfrac{\mu}{\sigma},\\
\tfrac{1}{\sigma} + 1\bigr)}, \quad k = 0, 1, 2, \ldots\$\$

The Waring is the beta-geometric: a geometric distribution whose success
probability carries a beta prior. It is the two-parameter long-tailed
count law behind accident proneness and repeat-buying models, and
generalises the
[Yule-Simon](https://janolefi.github.io/RTMBdist/reference/yules.md)
distribution, which is the case \\\sigma = \mu\\ shifted to start at
one. The variance \$\$\mathrm{Var}(X) = \frac{\mu (\sigma + 1)(\mu +
1)}{1 - \sigma}\$\$ is finite only for \\\sigma \< 1\\, and the tail is
a power law throughout.

It is exactly the [mean-parameterised beta-negative
binomial](https://janolefi.github.io/RTMBdist/reference/bnbinom2.md)
with `nu = 1`. Fixing `size` at one is what gives it a closed-form
distribution function, which the beta-negative binomial does not have in
general, so one-step-ahead residuals are available here but not there.

## References

Irwin, J. O. (1963) The place of mathematics in medical and biological
statistics. Journal of the Royal Statistical Society A, 126, 1-45,
doi:10.2307/2982445.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[yules](https://janolefi.github.io/RTMBdist/reference/yules.md),
[bnbinom2](https://janolefi.github.io/RTMBdist/reference/bnbinom2.md),
[bnbinom](https://janolefi.github.io/RTMBdist/reference/bnbinom.md)

## Examples

``` r
set.seed(123)
x <- rwaring(5, mu = 2, sigma = 0.5)
d <- dwaring(x, mu = 2, sigma = 0.5)
p <- pwaring(x, mu = 2, sigma = 0.5)
```
