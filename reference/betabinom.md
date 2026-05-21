# Beta-binomial distribution

Density and random generation for the beta-binomial distribution.

## Usage

``` r
dbetabinom(x, size, shape1, shape2, log = FALSE)

rbetabinom(n, size, shape1, shape2)
```

## Arguments

- x:

  vector of non-negative counts.

- size:

  vector of total counts (number of trials). Needs to be \>= `x`.

- shape1:

  positive shape parameter 1 of the Beta prior.

- shape2:

  positive shape parameter 2 of the Beta prior.

- log:

  logical; if `TRUE`, densities are returned on the log scale.

- n:

  number of random values to return (for `rbetabinom`).

## Value

`dbetabinom` gives the density and `rbetabinom` generates random
samples.

## Details

This implementation of `dbetabinom` allows for automatic differentiation
with `RTMB`.

\$\$P(X = k;\\ n, a, b) = \binom{n}{k} \frac{B(k+a,\\ n-k+b)}{B(a,\\
b)}, \quad k = 0, 1, \ldots, n.\$\$

## Examples

``` r
set.seed(123)
x <- rbetabinom(1, 10, 2, 5)
d <- dbetabinom(x, 10, 2, 5)
```
