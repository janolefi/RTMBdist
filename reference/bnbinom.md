# Beta-negative binomial distribution

Probability mass function and random generation for the beta-negative
binomial distribution.

## Usage

``` r
dbnbinom(x, size, shape1, shape2, log = FALSE)

rbnbinom(n, size, shape1, shape2)
```

## Arguments

- x:

  vector of non-negative counts.

- size:

  positive number of successes (need not be an integer).

- shape1:

  positive shape parameter 1 of the beta prior.

- shape2:

  positive shape parameter 2 of the beta prior.

- log:

  logical; if `TRUE`, probabilities are returned on the log scale.

- n:

  number of random values to return (for `rbnbinom`).

## Value

`dbnbinom` gives the probability mass function and `rbnbinom` generates
random deviates.

## Details

`dbnbinom` allows for automatic differentiation with `RTMB`.

The beta-negative binomial arises by giving the success probability of a
negative binomial a beta prior, in the same way that the
[beta-binomial](https://janolefi.github.io/RTMBdist/reference/betabinom.md)
does for the binomial: \$\$P(X = k;\\ r, a, b) = \frac{\Gamma(k +
r)}{k!\\ \Gamma(r)} \frac{B(a + r,\\ b + k)}{B(a,\\ b)}, \quad k = 0, 1,
2, \ldots\$\$

The extra beta layer gives a much heavier tail than the negative
binomial: the mean \\rb / (a - 1)\\ exists only for \\a \> 1\\ and the
variance \$\$\frac{r b (r + a - 1)(b + a - 1)}{(a - 2)(a - 1)^2}\$\$
only for \\a \> 2\\. As \\a \to \infty\\ with \\b/(a+b)\\ held fixed the
distribution collapses to the negative binomial.

Note that the three parameters are only weakly identified: the
likelihood has a long ridge along which `size` and `shape2` trade off
against each other, so fitting all three at once needs either a lot of
data or a restriction. The [mean
parameterisation](https://janolefi.github.io/RTMBdist/reference/bnbinom2.md)
is usually the more stable one to estimate in.

There is no distribution function, since the beta-negative binomial
distribution function has no closed form and the support is unbounded.
One-step-ahead residuals are therefore not available.

## References

Johnson, N. L., Kemp, A. W. and Kotz, S. (2005) Univariate Discrete
Distributions, 3rd edition, Wiley, doi:10.1002/0471715816.

## See also

[bnbinom2](https://janolefi.github.io/RTMBdist/reference/bnbinom2.md),
[betabinom](https://janolefi.github.io/RTMBdist/reference/betabinom.md),
[nbinom2](https://janolefi.github.io/RTMBdist/reference/nbinom2.md)

## Examples

``` r
set.seed(123)
x <- rbnbinom(5, size = 3, shape1 = 4, shape2 = 2)
d <- dbnbinom(x, size = 3, shape1 = 4, shape2 = 2)
```
