# Gaussian copula constructor

Returns a function computing the log density of the bivariate Gaussian
copula, intended to be used with
[`dcopula`](https://janolefi.github.io/RTMBdist/reference/dcopula.md).

## Usage

``` r
cgaussian(rho = 0)
```

## Arguments

- rho:

  Correlation parameter (\\-1 \< rho \< 1\\).

## Value

Function of two arguments (u,v) returning log copula density.

The Gaussian copula density is \$\$ c(u,v;\rho) =
\frac{1}{\sqrt{1-\rho^2}} \exp\left\\-\frac{1}{2(1-\rho^2)} (z_1^2 - 2
\rho z_1 z_2 + z_2^2) + \frac{1}{2}(z_1^2 + z_2^2) \right\\, \$\$ where
\\z_1 = \Phi^{-1}(u)\\, \\z_2 = \Phi^{-1}(v)\\, and \\-1 \< \rho \< 1\\.

## See also

[`cclayton()`](https://janolefi.github.io/RTMBdist/reference/cclayton.md),
[`cgumbel()`](https://janolefi.github.io/RTMBdist/reference/cgumbel.md),
[`cfrank()`](https://janolefi.github.io/RTMBdist/reference/cfrank.md)

## Examples

``` r
x <- c(0.5, 1); y <- c(1, 2)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
dcopula(d1, d2, p1, p2, copula = cgaussian(0.5), log = TRUE)
#> [1] -2.818014 -4.809862
```
