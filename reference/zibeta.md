# Zero-inflated beta distribution

Density, distribution function, and random generation for the
zero-inflated beta distribution.

## Usage

``` r
dzibeta(x, shape1, shape2, zeroprob = 0, log = FALSE)

pzibeta(q, shape1, shape2, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rzibeta(n, shape1, shape2, zeroprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- shape1, shape2:

  non-negative shape parameters of the beta distribution

- zeroprob:

  zero-inflation probability between 0 and 1.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dzibeta` gives the density, `pzibeta` gives the distribution function,
and `rzibeta` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

\$\$f(x;\\a,b,p_0) = p_0\\\mathbf{1}\[x=0\] +
(1-p_0)\\f\_{\mathrm{Beta}}(x;\\a,b)\\\mathbf{1}\[x\in(0,1)\],\$\$ where
\\p_0\\ is `zeroprob`.

## Examples

``` r
set.seed(123)
x <- rzibeta(1, 2, 2, 0.5)
d <- dzibeta(x, 2, 2, 0.5)
p <- pzibeta(x, 2, 2, 0.5)
```
