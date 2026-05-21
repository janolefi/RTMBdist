# Zero- and one-inflated beta distribution

Density, distribution function, and random generation for the
zero-one-inflated beta distribution.

## Usage

``` r
dzoibeta(x, shape1, shape2, zeroprob = 0, oneprob = 0, log = FALSE)

pzoibeta(q, shape1, shape2, zeroprob = 0, oneprob = 0,
         lower.tail = TRUE, log.p = FALSE)

rzoibeta(n, shape1, shape2, zeroprob = 0, oneprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- shape1, shape2:

  non-negative shape parameters of the beta distribution

- zeroprob:

  zero-inflation probability between 0 and 1.

- oneprob:

  one-inflation probability between 0 and 1.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dzoibeta` gives the density, `pzoibeta` gives the distribution
function, and `rzoibeta` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

\$\$f(x;\\a,b,p_0,p_1) = p_0\\\mathbf{1}\[x=0\] +
(1-p_0-p_1)\\f\_{\mathrm{Beta}}(x;\\a,b)\\\mathbf{1}\[x\in(0,1)\] +
p_1\\\mathbf{1}\[x=1\],\$\$ where \\p_0\\ = `zeroprob` and \\p_1\\ =
`oneprob`.

## Examples

``` r
set.seed(123)
x <- rzoibeta(1, 2, 2, 0.2, 0.3)
d <- dzoibeta(x, 2, 2, 0.2, 0.3)
p <- pzoibeta(x, 2, 2, 0.2, 0.3)
```
