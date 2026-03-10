# Clayton copula constructors

Construct a function that computes the log density or CDF of the
bivariate Clayton copula, intended to be used with
[`dcopula`](https://janolefi.github.io/RTMBdist/reference/dcopula.md).

## Usage

``` r
cclayton(theta)

Cclayton(theta)
```

## Arguments

- theta:

  Positive dependence parameter (\\\theta \> 0\\).

## Value

A function of two arguments (u,v) returning log copula **density**
(`cclayton`) or copula **CDF** (`Cclayton`).

## Details

The Clayton copula density is \$\$ c(u,v;\theta) = (1+\theta)
(uv)^{-(1+\theta)} \left( u^{-\theta} + v^{-\theta} - 1
\right)^{-(2\theta+1)/\theta}, \quad \theta \> 0. \$\$

## See also

[`cgaussian()`](https://janolefi.github.io/RTMBdist/reference/cgaussian.md),
[`cgumbel()`](https://janolefi.github.io/RTMBdist/reference/cgumbel.md),
[`cfrank()`](https://janolefi.github.io/RTMBdist/reference/cfrank.md)

## Examples

``` r
x <- c(0.5, 1); y <- c(0.2, 0.8)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
dcopula(d1, d2, p1, p2, copula = cclayton(2), log = TRUE)
#> [1] -3.8093660 -0.1671136

# CDF version (for discrete copulas)
Cclayton(1.5)(0.5, 0.4)
#> [1] 0.3104447
```
