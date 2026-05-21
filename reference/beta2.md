# Reparameterised beta distribution

Density, distribution function, quantile function, and random generation
for the beta distribution reparameterised in terms of mean and
concentration.

## Usage

``` r
dbeta(x, shape1, shape2, log = FALSE, eps = 0)

dbeta2(x, mu, phi, log = FALSE, eps = 0)

pbeta2(q, mu, phi, lower.tail = TRUE, log.p = FALSE)

qbeta2(p, mu, phi, lower.tail = TRUE, log.p = FALSE)

rbeta2(n, mu, phi)
```

## Arguments

- x, q:

  vector of quantiles

- shape1, shape2:

  non-negative parameters

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- eps:

  for internal use only, don't change.

- mu:

  mean parameter, must be in the interval from 0 to 1.

- phi:

  concentration parameter, must be positive.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return.

## Value

`dbeta2` gives the density, `pbeta2` gives the distribution function,
`qbeta2` gives the quantile function, and `rbeta2` generates random
deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

Currently, `dbeta` masks `RTMB::dbeta` because the latter has a
numerically unstable gradient.

`dbeta` and `dbeta2` compute the beta density. `dbeta2` reparameterises
by mean \\\mu\\ and concentration \\\phi\\: \$\$f(x; \mu, \phi) =
\frac{x^{\mu\phi - 1}(1-x)^{(1-\mu)\phi - 1}}{B(\mu\phi,\\
(1-\mu)\phi)}, \quad x \in (0,1).\$\$

## Examples

``` r
set.seed(123)
x <- rbeta2(1, 0.5, 1)
d <- dbeta2(x, 0.5, 1)
p <- pbeta2(x, 0.5, 1)
q <- qbeta2(p, 0.5, 1)
```
