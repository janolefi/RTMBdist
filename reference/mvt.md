# Multivariate t distribution

Density and and random generation for the multivariate t distribution

## Usage

``` r
dmvt(x, mu, Sigma, df, log = FALSE)

rmvt(n, mu, Sigma, df)
```

## Arguments

- x:

  vector or matrix of quantiles

- mu:

  vector or matrix of location parameters (mean if `df` \> 1)

- Sigma:

  positive definite scale matrix (proportional to the covariance matrix
  if `df` \> 2)

- df:

  degrees of freedom; must be positive

- log:

  logical; if `TRUE`, densities \\p\\ are returned as \\\log(p)\\.

- n:

  number of random values to return.

## Value

`dmvt` gives the density, `rmvt` generates random deviates.

## Details

This implementation of `dmvt` allows for automatic differentiation with
`RTMB`.

Note: for `df` \\\le 1\\ the mean is undefined, and for `df` \\\le 2\\
the covariance is infinite. For `df` \> 2, the covariance is
`df/(df-2) * Sigma`.

\$\$f(\mathbf{x};\\\boldsymbol{\mu},\Sigma,\nu) =
\frac{\Gamma((\nu+d)/2)}{\Gamma(\nu/2)\\(\nu\pi)^{d/2}\\\|\Sigma\|^{1/2}}
\left(1 + \frac{(\mathbf{x}-\boldsymbol{\mu})^\top
\Sigma^{-1}(\mathbf{x}-\boldsymbol{\mu})}{\nu}\right)^{-(\nu+d)/2},\$\$
where \\d\\ is the dimension.

## Examples

``` r
# single mu
mu <- c(1,2,3)
Sigma <- diag(c(1,1,1))
df <- 5
x <- rmvt(2, mu, Sigma, df)
d <- dmvt(x, mu, Sigma, df)
# vectorised over mu
mu <- rbind(c(1,2,3), c(0, 0.5, 1))
x <- rmvt(2, mu, Sigma, df)
d <- dmvt(x, mu, Sigma, df)
```
