# Zero-inflated and reparameterised gamma distribution

Density, distribution function, and random generation for the
zero-inflated gamma distribution reparameterised in terms of mean and
standard deviation.

## Usage

``` r
dzigamma2(x, mean = 1, sd = 1, zeroprob = 0, log = FALSE)

pzigamma2(q, mean = 1, sd = 1, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rzigamma2(n, mean = 1, sd = 1, zeroprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- mean:

  mean parameter, must be positive.

- sd:

  standard deviation parameter, must be positive.

- zeroprob:

  zero-inflation probability between 0 and 1.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- n:

  number of random values to return

## Value

`dzigamma2` gives the density, `pzigamma2` gives the distribution
function, and `rzigamma` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

## Examples

``` r
x <- rzigamma2(1, 2, 1, 0.5)
d <- dzigamma2(x, 2, 1, 0.5)
p <- pzigamma2(x, 2, 1, 0.5)
```
