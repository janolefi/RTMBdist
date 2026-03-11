# Beta prime distribution

Density, distribution function, quantile function, and random generation
for the Beta prime distribution.

## Usage

``` r
dbetaprime(x, shape1, shape2, log = FALSE)

pbetaprime(q, shape1, shape2, lower.tail = TRUE, log.p = FALSE)

qbetaprime(p, shape1, shape2, lower.tail = TRUE, log.p = FALSE)

rbetaprime(n, shape1, shape2)
```

## Arguments

- x, q:

  vector of quantiles

- shape1, shape2:

  non-negative shape parameters of the corresponding Beta distribution

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return.

## Value

`dbetaprime` gives the density, `pbetaprime` gives the distribution
function, `qbetaprime` gives the quantile function, and `rbetaprime`
generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

If \\X \sim \text{Beta}(\alpha, \beta)\\, then \\\frac{X}{1-X} \sim
\text{Betaprime}(\alpha, \beta)\\

## Examples

``` r
x <- rbetaprime(1, 2, 1)
d <- dbetaprime(x, 2, 1)
p <- pbetaprime(x, 2, 1)
q <- qbetaprime(p, 2, 1)
```
