# Gumbel distribution

Density, distribution function, quantile function, and random generation
for the Gumbel distribution.

## Usage

``` r
dgumbel(x, location = 0, scale = 1, log = FALSE)

pgumbel(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)

qgumbel(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)

rgumbel(n, location = 0, scale = 1)
```

## Arguments

- x, q:

  vector of quantiles

- location:

  location parameter

- scale:

  scale parameter, must be positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dgumbel` gives the density, `pgumbel` gives the distribution function,
`qgumbel` gives the quantile function, and `rgumbel` generates random
deviates.

## Details

This implementation of `dgumbel` allows for automatic differentiation
with `RTMB`.

\$\$f(x;\\\mu,\sigma) = \frac{1}{\sigma}\exp\\\bigl(-(z + e^{-z})\bigr),
\quad z = \frac{x-\mu}{\sigma}.\$\$

## Examples

``` r
x <- rgumbel(1, 0.5, 2)
d <- dgumbel(x, 0.5, 2)
p <- pgumbel(x, 0.5, 2)
q <- qgumbel(p, 0.5, 2)
```
