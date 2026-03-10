# Moment-parameterised skew t distribution

Density, distribution function, quantile function, and random generation
for the skew t distribution reparameterised so that `mean` and `sd`
correspond to the *true* mean and standard deviation.

## Usage

``` r
dskewt2(x, mean = 0, sd = 1, skew = 0, df = 100, log = FALSE)

pskewt2(q, mean = 0, sd = 1, skew = 0, df = 100,
        method = 0, lower.tail = TRUE, log.p = FALSE)

qskewt2(p, mean = 0, sd = 1, skew = 0, df = 100, tol = 1e-08, method = 0)

rskewt2(n, mean = 0, sd = 1, skew = 0, df = 100)
```

## Arguments

- x, q:

  vector of quantiles

- mean:

  mean parameter

- sd:

  standard deviation parameter, must be positive.

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

`dskewt2` gives the density, `pskewt2` gives the distribution function,
`qskewt2` gives the quantile function, and `rskewt2` generates random
deviates.

## Details

This corresponds to the skew t type 2 distribution in GAMLSS
([`ST2`](https://rdrr.io/pkg/gamlss.dist/man/ST1.html)), see pp. 411-412
of Rigby et al. (2019) and the version implemented in the `sn` package.
However, it is reparameterised in terms of a standard deviation
parameter `sd` rather than just a scale parameter `sigma`. Details of
this reparameterisation are given below. This implementation of `dskewt`
allows for automatic differentiation with `RTMB` while the other
functions are imported from the `sn` package. See
`sn::`[`dst`](https://rdrr.io/pkg/sn/man/dst.html) for more details.

**Caution:** In a numerial optimisation, the `skew` parameter should
NEVER be initialised with exactly zero. This will cause the initial and
all subsequent derivatives to be exactly zero and hence the parameter
will remain at its initial value.

For given `skew` \\= \alpha\\ and `df` = \\\nu\\, define \$\$ \delta =
\alpha / \sqrt{1 + \alpha^2}, \qquad b\_\nu = \sqrt{\nu / \pi}\\
\Gamma((\nu-1)/2)/\Gamma(\nu/2), \$\$ then \$\$ E(X) = \mu + \sigma
\delta b\_\nu,\quad Var(X) = \sigma^2 \left( \frac{\nu}{\nu-2} -
\delta^2 b\_\nu^2 \right). \$\$

## See also

[skewt](https://janolefi.github.io/RTMBdist/reference/skewt.md),
[skewnorm](https://janolefi.github.io/RTMBdist/reference/skewnorm.md),
[skewnorm2](https://janolefi.github.io/RTMBdist/reference/skewnorm2.md)

## Examples

``` r
x <- rskewt2(1, 1, 2, 5, 5)
d <- dskewt2(x, 1, 2, 5, 5)
p <- pskewt2(x, 1, 2, 5, 5)
q <- qskewt2(p, 1, 2, 5, 5)
```
