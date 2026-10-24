# Birnbaum-Saunders distribution

Density, distribution function, quantile function, and random generation
for the Birnbaum-Saunders (fatigue life) distribution.

## Usage

``` r
dbisa(x, scale = 1, shape, log = FALSE)

pbisa(q, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)

qbisa(p, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)

rbisa(n, scale = 1, shape)
```

## Arguments

- x, q:

  vector of quantiles.

- scale:

  scale parameter \\\beta\\, the median, must be positive.

- shape:

  shape parameter \\\alpha\\, must be positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities.

- n:

  number of random values to return.

## Value

`dbisa` gives the density, `pbisa` gives the distribution function,
`qbisa` gives the quantile function, and `rbisa` generates random
deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

The Birnbaum-Saunders distribution with scale \\\beta\\ and shape
\\\alpha\\ is the distribution of \\X\\ for which \$\$\xi(X) =
\frac{1}{\alpha}\left(\sqrt{X/\beta} - \sqrt{\beta/X}\right) =
\frac{2}{\alpha}\sinh\left(\tfrac{1}{2}\log(X/\beta)\right)\$\$ is
standard normal, so that \\F(x) = \Phi(\xi(x))\\ and \$\$f(x) =
\phi(\xi(x))\\ \frac{x + \beta}{2\alpha\sqrt{\beta}\\ x^{3/2}}, \quad x
\> 0.\$\$ The scale \\\beta\\ is the median, and \\E(X) = \beta(1 +
\alpha^2/2)\\. The quantile function is \\\beta
\exp(2\\\mathrm{asinh}(\alpha z / 2))\\ with \\z = \Phi^{-1}(p)\\.

## References

Birnbaum, Z. W. and Saunders, S. C. (1969) A new family of life
distributions. Journal of Applied Probability, 6, 319-327.

## See also

[invgauss](https://janolefi.github.io/RTMBdist/reference/invgauss.md),
[pgweibull](https://janolefi.github.io/RTMBdist/reference/pgweibull.md),
[gompertz](https://janolefi.github.io/RTMBdist/reference/gompertz.md)

## Examples

``` r
set.seed(1)
x <- rbisa(5, scale = 2, shape = 0.5)
d <- dbisa(x, 2, 0.5)
p <- pbisa(x, 2, 0.5)
q <- qbisa(p, 2, 0.5)
```
