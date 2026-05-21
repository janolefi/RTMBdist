# Reparameterised zero- and one-inflated beta distribution

Density, distribution function, and random generation for the
zero-one-inflated beta distribution reparameterised in terms of mean and
concentration.

## Usage

``` r
dzoibeta2(x, mu, phi, zeroprob = 0, oneprob = 0, log = FALSE)

pzoibeta2(q, mu, phi, zeroprob = 0, oneprob = 0,
         lower.tail = TRUE, log.p = FALSE)

rzoibeta2(n, mu, phi, zeroprob = 0, oneprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  mean parameter, must be in the interval from 0 to 1.

- phi:

  concentration parameter, must be positive.

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

`dzoibeta2` gives the density, `pzoibeta2` gives the distribution
function, and `rzoibeta2` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

Uses the same density as `zoibeta` with \\a = \mu\phi\\ and \\b =
(1-\mu)\phi\\.

## Examples

``` r
set.seed(123)
x <- rzoibeta2(1, 0.6, 2, 0.2, 0.3)
d <- dzoibeta2(x, 0.6, 2, 0.2, 0.3)
p <- pzoibeta2(x, 0.6, 2, 0.2, 0.3)
```
