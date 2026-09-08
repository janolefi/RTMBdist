# Zero-truncated geometric distribution

Probability mass function, distribution function, and random generation
for the zero-truncated geometric distribution.

## Usage

``` r
dztgeom(x, prob, log = FALSE)

pztgeom(q, prob, lower.tail = TRUE, log.p = FALSE)

rztgeom(n, prob)
```

## Arguments

- x, q:

  integer vector of counts

- prob:

  probability of success in each trial, in (0,1)

- log, log.p:

  logical; return log-density if TRUE

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dztgeom` gives the probability mass function, `pztgeom` gives the
distribution function, and `rztgeom` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

By definition, this distribution only has support on the positive
integers (1, 2, ...). Any zero-truncated distribution is defined as
\$\$P(X=x \| X\>0) = P(X=x) / (1 - P(X=0)),\$\$ where \\P(X=x)\\ is the
probability mass function of the corresponding untruncated distribution.
For the geometric with success probability \\\pi\\ this gives \$\$P(X=x
\| X\>0) = \pi\\(1-\pi)^{x-1}, \quad x = 1, 2, \ldots\$\$

## See also

[zigeom](https://janolefi.github.io/RTMBdist/reference/zigeom.md),
[hgeom](https://janolefi.github.io/RTMBdist/reference/hgeom.md),
[ztnbinom](https://janolefi.github.io/RTMBdist/reference/ztnbinom.md),
[ztpois](https://janolefi.github.io/RTMBdist/reference/ztpois.md)

## Examples

``` r
set.seed(123)
x <- rztgeom(5, prob = 0.3)
d <- dztgeom(x, prob = 0.3)
p <- pztgeom(x, prob = 0.3)
```
