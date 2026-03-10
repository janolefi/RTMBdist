# Frank copula constructor

Returns a function computing the log density of the bivariate Frank
copula, intended to be used with
[`dcopula`](https://janolefi.github.io/RTMBdist/reference/dcopula.md).

## Usage

``` r
cfrank(theta)

Cfrank(theta)
```

## Arguments

- theta:

  Dependence parameter (\\\theta \neq 0\\).

## Value

A function of two arguments `(u, v)` returning either the log copula
density (`cfrank`) or the copula CDF (`Cfrank`).

## Details

The Frank copula density is \$\$ c(u,v;\theta) = \frac{\theta
(1-e^{-\theta}) e^{-\theta(u+v)}} {\left\[(e^{-\theta u}-1)(e^{-\theta
v}-1) + (1 - e^{-\theta}) \right\]^2}, \quad \theta \ne 0. \$\$

## See also

[`cgaussian()`](https://janolefi.github.io/RTMBdist/reference/cgaussian.md),
[`cclayton()`](https://janolefi.github.io/RTMBdist/reference/cclayton.md),
[`cgumbel()`](https://janolefi.github.io/RTMBdist/reference/cgumbel.md)

## Examples

``` r
x <- c(0.5, 1); y <- c(1, 2)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
dcopula(d1, d2, p1, p2, copula = cfrank(2), log = TRUE)
#> [1] -4.585248 -7.325831
```
