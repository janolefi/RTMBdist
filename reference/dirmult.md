# Dirichlet-multinomial distribution

Density and and random generation for the Dirichlet-multinomial
distribution.

## Usage

``` r
ddirmult(x, size, alpha, log = FALSE)

rdirmult(n, size, alpha)
```

## Arguments

- x:

  vector or matrix of non-negative counts, where rows are observations
  and columns are categories.

- size:

  vector of total counts for each observation. Needs to match the row
  sums of `x`.

- alpha:

  vector or matrix of positive shape parameters

- log:

  logical; if `TRUE`, densities \\p\\ are returned as \\\log(p)\\.

- n:

  number of random values to return.

## Value

`ddirmult` gives the density and `rdirmult` generates random samples.

## Details

This implementation of `ddirmult` allows for automatic differentiation
with `RTMB`.

\$\$P(\mathbf{x};\\\boldsymbol{\alpha}, n) =
\frac{\Gamma(\alpha_0)\\\Gamma(n+1)}{\Gamma(n+\alpha_0)} \prod_i
\frac{\Gamma(x_i + \alpha_i)}{\Gamma(\alpha_i)\\\Gamma(x_i+1)},\$\$
where \\\alpha_0 = \sum_i \alpha_i\\ and \\n = \sum_i x_i\\.

## Examples

``` r
# single alpha
alpha <- c(1,2,3)
size <- 10
x <- rdirmult(1, size, alpha)
d <- ddirmult(x, size, alpha)
# vectorised over alpha and size
alpha <- rbind(alpha, 2*alpha)
size <- c(size, 3*size)
x <- rdirmult(2, size, alpha)
```
