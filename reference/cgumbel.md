# Gumbel copula constructors

Construct functions that compute either the log density or the CDF of
the bivariate Gumbel copula, intended for use with
[`dcopula`](https://janolefi.github.io/RTMBdist/reference/dcopula.md).

## Usage

``` r
cgumbel(theta)

Cgumbel(theta)
```

## Arguments

- theta:

  Dependence parameter (\\\theta \>= 1\\).

## Value

A function of two arguments `(u, v)` returning either the log copula
density (`cgumbel`) or the copula CDF (`Cgumbel`).

## Details

The Gumbel copula density

\$\$ c(u,v;\theta) = \frac{e^{-A}}{uv} \\ (\tilde u \tilde v)^{\theta -
1} S^{1/\theta - 2} (A + \theta - 1), \$\$ where \\\tilde u = -\log u\\,
\\\tilde v = -\log v\\, \\S = \tilde u^\theta + \tilde v^\theta\\ and
\\A = S^{1/\theta}\\, so that \\e^{-A}\\ is the copula CDF.

## See also

[`cgaussian()`](https://janolefi.github.io/RTMBdist/reference/cgaussian.md),
[`cclayton()`](https://janolefi.github.io/RTMBdist/reference/cclayton.md),
[`cfrank()`](https://janolefi.github.io/RTMBdist/reference/cfrank.md)

## Examples

``` r
x <- c(0.5, 1); y <- c(0.2, 0.4)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
dcopula(d1, d2, p1, p2, copula = cgumbel(1.5), log = TRUE)
#> [1] -1.751222 -1.135706

# CDF version (for discrete copulas)
Cgumbel(1.5)(0.5, 0.4)
#> [1] 0.2770518
```
