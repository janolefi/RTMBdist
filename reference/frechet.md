# Frechet distribution

Density, distribution function, quantile function, and random generation
for the Frechet distribution.

## Usage

``` r
dfrechet(x, mu = 0, sigma = 1, alpha = 1, log = FALSE)

pfrechet(q, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE)

qfrechet(p, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE)

rfrechet(n, mu = 0, sigma = 1, alpha = 1)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter, the lower end point of the support.

- sigma:

  scale parameter, must be positive.

- alpha:

  shape parameter, must be positive.

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

`dfrechet` gives the density, `pfrechet` gives the distribution
function, `qfrechet` gives the quantile function, and `rfrechet`
generates random deviates.

## Details

`dfrechet` and `pfrechet` allow for automatic differentiation with
`RTMB`.

With \\y = (x - \mu) / \sigma\\ the density is
\$\$f(x;\\\mu,\sigma,\alpha) = \frac{\alpha}{\sigma} y^{-\alpha - 1}
\exp(-y^{-\alpha}), \quad x \> \mu,\$\$ and the distribution function is
\\F(x) = \exp(-y^{-\alpha})\\.

The Frechet distribution is the heavy-tailed extreme value distribution:
it is the [generalised extreme
value](https://janolefi.github.io/RTMBdist/reference/gev.md)
distribution with shape \\\xi = 1/\alpha \> 0\\, reparameterised so that
the shape enters as a tail index rather than as a reciprocal. All
moments of order \\\alpha\\ and above are infinite, so the mean exists
only for \\\alpha \> 1\\ and the variance only for \\\alpha \> 2\\.

## References

Frechet, M. (1927) Sur la loi de probabilite de l'ecart maximum. Annales
de la Societe Polonaise de Mathematique, 6, 93-116.

Kotz, S. and Nadarajah, S. (2000) Extreme Value Distributions: Theory
and Applications, Imperial College Press, doi:10.1142/p191.

## See also

[gev](https://janolefi.github.io/RTMBdist/reference/gev.md),
[gumbel](https://janolefi.github.io/RTMBdist/reference/gumbel.md),
[pareto](https://janolefi.github.io/RTMBdist/reference/pareto.md)

## Examples

``` r
set.seed(123)
x <- rfrechet(5, mu = 0, sigma = 1, alpha = 3)
d <- dfrechet(x, mu = 0, sigma = 1, alpha = 3)
p <- pfrechet(x, mu = 0, sigma = 1, alpha = 3)
q <- qfrechet(p, mu = 0, sigma = 1, alpha = 3)
```
