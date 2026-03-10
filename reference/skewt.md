# Skewed students t distribution

Density, distribution function, quantile function, and random generation
for the skew t distribution (type 2).

## Usage

``` r
dskewt(x, mu = 0, sigma = 1, skew = 0, df = 100, log = FALSE)

pskewt(q, mu = 0, sigma = 1, skew = 0, df = 100,
       method = 0, lower.tail = TRUE, log.p = FALSE)

qskewt(p, mu = 0, sigma = 1, skew = 0, df = 100,
       tol = 1e-8, method = 0)

rskewt(n, mu = 0, sigma = 1, skew = 0, df = 100)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter

- sigma:

  scale parameter, must be positive.

- skew:

  skewness parameter, can be positive or negative.

- df:

  degrees of freedom, must be positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- method:

  an integer value between 0 and 5 which selects the computing method;
  see ‘Details’ in the [`pst`](https://rdrr.io/pkg/sn/man/dst.html)
  documentation below for the meaning of these values. If method=0
  (default value), an automatic choice is made among the four actual
  computing methods, depending on the other arguments.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- tol:

  a scalar value which regulates the accuracy of the result of qsn,
  measured on the probability scale.

- n:

  number of random values to return.

## Value

`dskewt` gives the density, `pskewt` gives the distribution function,
`qskewt` gives the quantile function, and `rskewt` generates random
deviates.

## Details

This corresponds to the skew t type 2 distribution in GAMLSS
([`ST2`](https://rdrr.io/pkg/gamlss.dist/man/ST1.html)), see pp. 411-412
of Rigby et al. (2019) and the version implemented in the `sn` package.
This implementation of `dskewt` allows for automatic differentiation
with `RTMB` while the other functions are imported from the `sn`
package. See `sn::`[`dst`](https://rdrr.io/pkg/sn/man/dst.html) for more
details.

**Caution:** In a numerial optimisation, the `skew` parameter should
NEVER be initialised with exactly zero. This will cause the initial and
all subsequent derivatives to be exactly zero and hence the parameter
will remain at its initial value.

## See also

[skewt2](https://janolefi.github.io/RTMBdist/reference/skewt2.md),
[skewnorm](https://janolefi.github.io/RTMBdist/reference/skewnorm.md),
[skewnorm2](https://janolefi.github.io/RTMBdist/reference/skewnorm2.md)

## Examples

``` r
x <- rskewt(1, 1, 2, 5, 2)
d <- dskewt(x, 1, 2, 5, 2)
p <- pskewt(x, 1, 2, 5, 2)
q <- qskewt(p, 1, 2, 5, 2)
```
