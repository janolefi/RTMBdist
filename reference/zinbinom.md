# Zero-inflated negative binomial distribution

Probability mass function, distribution function, quantile function, and
random generation for the zero-inflated negative binomial distribution.

## Usage

``` r
dzinbinom(x, size, prob, zeroprob = 0, log = FALSE)

pzinbinom(q, size, prob, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rzinbinom(n, size, prob, zeroprob = 0)
```

## Arguments

- x, q:

  vector of (non-negative integer) quantiles

- size:

  size parameter, must be positive.

- prob:

  mean parameter, must be positive.

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

- p:

  vector of probabilities

## Value

`dzinbinom` gives the density, `pzinbinom` gives the distribution
function, and `rzinbinom` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

\$\$P(X=k;\\r,p,p_0) = p_0\\\mathbf{1}\[k=0\] +
(1-p_0)\\P\_{\mathrm{NB}}(k;\\r,p),\$\$ where \\p_0\\ is `zeroprob`.

## Examples

``` r
set.seed(123)
x <- rzinbinom(1, size = 2, prob = 0.5, zeroprob = 0.5)
d <- dzinbinom(x, size = 2, prob = 0.5, zeroprob = 0.5)
p <- pzinbinom(x, size = 2, prob = 0.5, zeroprob = 0.5)
```
