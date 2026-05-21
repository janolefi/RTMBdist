# Zero-inflated gamma distribution

Density, distribution function, and random generation for the
zero-inflated gamma distribution.

## Usage

``` r
dzigamma(x, shape, scale, zeroprob = 0, log = FALSE)

pzigamma(q, shape, scale, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rzigamma(n, shape, scale, zeroprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- shape:

  positive shape parameter

- scale:

  positive scale parameter

- zeroprob:

  zero-inflation probability between 0 and 1.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return

## Value

`dzigamma` gives the density, `pzigamma` gives the distribution
function, and `rzigamma` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

\$\$f(x;\\\alpha,s,p_0) = p_0\\\mathbf{1}\[x=0\] +
(1-p_0)\\f\_{\mathrm{Gamma}}(x;\\\alpha,s)\\\mathbf{1}\[x\>0\],\$\$
where \\p_0\\ is `zeroprob`.

## Examples

``` r
x <- rzigamma(1, 1, 1, 0.5)
d <- dzigamma(x, 1, 1, 0.5)
p <- pzigamma(x, 1, 1, 0.5)
```
