# Folded normal distribution

Density, distribution function, and random generation for the folded
normal distribution.

## Usage

``` r
dfoldnorm(x, mu = 0, sigma = 1, log = FALSE)

pfoldnorm(q, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)

rfoldnorm(n, mu = 0, sigma = 1)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter

- sigma:

  scale parameter, must be positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return

- p:

  vector of probabilities

## Value

`dfoldnorm` gives the density, `pfoldnorm` gives the distribution
function, and `rfoldnorm` generates random deviates.

## Details

This implementation of `dfoldnorm` allows for automatic differentiation
with `RTMB`.

\$\$f(x;\\\mu,\sigma) =
\frac{1}{\sigma\sqrt{2\pi}}\left\[\exp\\\left(-\frac{(x-\mu)^2}{2\sigma^2}\right) +
\exp\\\left(-\frac{(x+\mu)^2}{2\sigma^2}\right)\right\], \quad x \geq
0.\$\$

With \\\mu = 0\\ this is the half-normal distribution, the limiting case
of the [half-t](https://janolefi.github.io/RTMBdist/reference/halft.md)
as its degrees of freedom grow.

## See also

[halft](https://janolefi.github.io/RTMBdist/reference/halft.md),
[halfcauchy](https://janolefi.github.io/RTMBdist/reference/halfcauchy.md)

## Examples

``` r
x <- rfoldnorm(1, 1, 2)
d <- dfoldnorm(x, 1, 2)
p <- pfoldnorm(x, 1, 2)
```
