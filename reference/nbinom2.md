# Reparameterised negative binomial distribution

Probability mass function, distribution function, quantile function, and
random generation for the negative binomial distribution reparameterised
in terms of mean and size.

## Usage

``` r
dnbinom2(x, mu, size, log = FALSE)

pnbinom2(q, mu, size, lower.tail = TRUE, log.p = FALSE)

qnbinom2(p, mu, size, lower.tail = TRUE, log.p = FALSE)

rnbinom2(n, mu, size)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  mean parameter, must be positive.

- size:

  size parameter, must be positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return.

- prob:

  probability of success in each trial. 0 \< prob \<= 1.

## Value

`dnbinom2` gives the density, `pnbinom2` gives the distribution
function, `qnbinom2` gives the quantile function, and `rnbinom2`
generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

`pnbinom` is an AD-compatible implementation of the standard
parameterisation of the CDF, missing from `RTMB`.

## Examples

``` r
set.seed(123)
x <- rnbinom2(1, 1, 2)
d <- dnbinom2(x, 1, 2)
p <- pnbinom2(x, 1, 2)
q <- qnbinom2(p, 1, 2)
```
