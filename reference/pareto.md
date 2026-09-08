# Pareto distribution

Density, distribution function, quantile function, and random generation
for the pareto distribution.

## Usage

``` r
dpareto(x, mu = 1, log = FALSE)

ppareto(q, mu = 1, lower.tail = TRUE, log.p = FALSE)

qpareto(p, mu = 1, lower.tail = TRUE, log.p = FALSE)

rpareto(n, mu = 1)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter, must be positive.

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

`dpareto` gives the density, `ppareto` gives the distribution function,
`qpareto` gives the quantile function, and `rpareto` generates random
deviates.

## Details

`dpareto` and `ppareto` allow for automatic differentiation with `RTMB`.
The parameterisation follows the `PARETO` family of the `gamlss.dist`
package.

\$\$f(x;\\\mu) = \frac{\mu}{x^{\mu+1}}, \quad x \> 1.\$\$

## References

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[gpd](https://janolefi.github.io/RTMBdist/reference/gpd.md),
[llogis](https://janolefi.github.io/RTMBdist/reference/llogis.md),
[gengamma](https://janolefi.github.io/RTMBdist/reference/gengamma.md)

## Examples

``` r
set.seed(123)
x <- rpareto(1, mu = 5)
d <- dpareto(x, mu = 5)
p <- ppareto(x, mu = 5)
q <- qpareto(p, mu = 5)
```
