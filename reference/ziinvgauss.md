# Zero-inflated inverse Gaussian distribution

Density, distribution function, and random generation for the
zero-inflated inverse Gaussian distribution.

## Usage

``` r
dziinvgauss(x, mean = 1, shape = 1, zeroprob = 0, log = FALSE)

pziinvgauss(q, mean = 1, shape = 1, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rziinvgauss(n, mean = 1, shape = 1, zeroprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- mean:

  location parameter

- shape:

  shape parameter, must be positive.

- zeroprob:

  zero-probability, must be in \\\[0, 1\]\\.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return

## Value

`dziinvgauss` gives the density, `pziinvgauss` gives the distribution
function, and `rziinvgauss` generates random deviates.

## Details

This implementation of `zidinvgauss` allows for automatic
differentiation with `RTMB`.

\$\$f(x;\\\mu,\lambda,p_0) = p_0\\\mathbf{1}\[x=0\] +
(1-p_0)\\f\_{\mathrm{IG}}(x;\\\mu,\lambda)\\\mathbf{1}\[x\>0\],\$\$
where \\p_0\\ is `zeroprob` and \\f\_{\mathrm{IG}}\\ is the inverse
Gaussian density.

## Examples

``` r
x <- rziinvgauss(1, 1, 2, 0.5)
d <- dziinvgauss(x, 1, 2, 0.5)
p <- pziinvgauss(x, 1, 2, 0.5)
```
