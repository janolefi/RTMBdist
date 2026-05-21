# Zero-inflated and reparameterised negative binomial distribution

Probability mass function, distribution function, quantile function and
random generation for the zero-inflated negative binomial distribution
reparameterised in terms of mean and size.

## Usage

``` r
dzinbinom2(x, mu, size, zeroprob = 0, log = FALSE)

pzinbinom2(q, mu, size, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rzinbinom2(n, mu, size, zeroprob = 0)
```

## Arguments

- x, q:

  vector of (non-negative integer) quantiles

- mu:

  mean parameter, must be positive.

- size:

  size parameter, must be positive.

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

`dzinbinom2` gives the density, `pzinbinom2` gives the distribution
function, and `rzinbinom2` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

Uses the same density as `zinbinom` with \\p = r/(r+\mu)\\:
\$\$P(X=k;\\\mu,r,p_0) = p_0\\\mathbf{1}\[k=0\] +
(1-p_0)\\P\_{\mathrm{NB}}\\\left(k;\\r,\tfrac{r}{r+\mu}\right).\$\$

## Examples

``` r
set.seed(123)
x <- rzinbinom2(1, 2, 1, zeroprob = 0.5)
d <- dzinbinom2(x, 2, 1, zeroprob = 0.5)
p <- pzinbinom2(x, 2, 1, zeroprob = 0.5)
```
