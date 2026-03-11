# Inverse Chi-squared distribution

Density, distribution function, quantile function, and random generation
for the inverse Chi-squared distribution.

## Usage

``` r
dinvchisq(x, df, scale = 1/df, log = FALSE)

pinvchisq(q, df, scale = 1/df, lower.tail = TRUE, log.p = FALSE)

qinvchisq(p, df, scale = 1/df, lower.tail = TRUE, log.p = FALSE)

rinvchisq(n, df, scale = 1/df)
```

## Arguments

- x, q:

  vector of quantiles, must be positive.

- df:

  degrees of freedom (\\\nu \> 0\\)

- scale:

  optional positive scale parameter. Default value of `1/df` corresponds
  to standard inverse gamma

- log, log.p:

  logical; if `TRUE`, probabilities/densities are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dinvchisq` gives the density, `pinvchisq` gives the distribution
function, `qinvchisq` gives the quantile function, and `rinvchisq`
generates random deviates.

## Details

If \\X \sim \text{Chisq}(\nu)\\, then \\1/X \sim \text{invChisq}(\nu)\\.

The inverse Chi-squared distribution with \\\nu\\ degrees of freedom has
density \$\$f(x) = \frac{(\nu/2)^{\nu/2}}{\Gamma(\nu/2)} x^{-(\nu/2+1)}
\exp(-\nu/(2x)), \quad x\>0.\$\$

This implementation of `dinvchisq`, `pinvchisq`, and `qinvchisq` allows
for automatic differentiation with `RTMB`.

## Examples

``` r
x <- rinvchisq(1, df = 5)
d <- dinvchisq(x, df = 5)
p <- pinvchisq(x, df = 5)
q <- qinvchisq(p, df = 5)
```
