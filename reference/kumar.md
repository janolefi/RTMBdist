# Kumaraswamy distribution

Density, distribution function, quantile function, and random generation
for the Kumaraswamy distribution.

## Usage

``` r
dkumar(x, a, b, log = FALSE)

pkumar(q, a, b, lower.tail = TRUE, log.p = FALSE)

qkumar(p, a, b, lower.tail = TRUE, log.p = FALSE)

rkumar(n, a, b)
```

## Arguments

- x, q:

  vector of quantiles in \\(0,1)\\

- a, b:

  positive shape parameters

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

`dkumar` gives the density, `pkumar` gives the distribution function,
`qkumar` gives the quantile function, and `rkumar` generates random
deviates.

## Examples

``` r
x <- rkumar(1, a = 1, b = 2)
d <- dkumar(x, a = 1, b = 2)
p <- pkumar(x, a = 1, b = 2)
q <- qkumar(p, a = 1, b = 2)
```
