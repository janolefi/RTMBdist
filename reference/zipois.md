# Zero-inflated Poisson distribution

Probability mass function, distribution function, and random generation
for the zero-inflated Poisson distribution.

## Usage

``` r
dzipois(x, lambda, zeroprob = 0, log = FALSE)

pzipois(q, lambda, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rzipois(n, lambda, zeroprob = 0)
```

## Arguments

- x, q:

  integer vector of counts

- lambda:

  vector of (non-negative) means

- zeroprob:

  zero-inflation probability between 0 and 1

- log, log.p:

  logical; return log-density if TRUE

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dzipois` gives the probability mass function, `pzipois` gives the
distribution function, and `rzipois` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

\$\$P(X=k;\\\lambda,p_0) = p_0\\\mathbf{1}\[k=0\] +
(1-p_0)\\P\_{\mathrm{Pois}}(k;\\\lambda),\$\$ where \\p_0\\ is
`zeroprob`.

## Examples

``` r
set.seed(123)
x <- rzipois(1, 0.5, 1)
d <- dzipois(x, 0.5, 1)
p <- pzipois(x, 0.5, 1)
```
