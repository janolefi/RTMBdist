# Joint density under a bivariate copula

Computes the joint density (or log-density) of a bivariate distribution
constructed from two arbitrary margins combined with a specified copula.

## Usage

``` r
dcopula(d1, d2, p1, p2, copula = cgaussian(0), log = FALSE)
```

## Arguments

- d1, d2:

  Marginal density values. If `log = TRUE`, supply the log-density. If
  `log = FALSE`, supply the raw density.

- p1, p2:

  Marginal CDF values. Need not be supplied on log scale.

- copula:

  A function of two arguments (`u`, `v`) returning the log copula
  density \\\log c(u,v)\\. You can either construct this yourself or use
  the copula constructors available (see details)

- log:

  Logical; if `TRUE`, return the log joint density. In this case, `d1`
  and `d2` must be on the log scale.

## Value

Joint density (or log-density) under the bivariate copula.

## Details

The joint density is \$\$f(x,y) = c(F_1(x), F_2(y)) \\ f_1(x)
f_2(y),\$\$ where \\F_i\\ are the marginal CDFs, \\f_i\\ are the
marginal densities, and \\c\\ is the copula density.

The marginal densities `d1`, `d2` and CDFs `p1`, `p2` must be
differentiable for automatic differentiation (AD) to work.

Available copula constructors are:

- [`cgaussian`](https://janolefi.github.io/RTMBdist/reference/cgaussian.md)
  (Gaussian copula)

- [`cclayton`](https://janolefi.github.io/RTMBdist/reference/cclayton.md)
  (Clayton copula)

- [`cgumbel`](https://janolefi.github.io/RTMBdist/reference/cgumbel.md)
  (Gumbel copula)

- [`cfrank`](https://janolefi.github.io/RTMBdist/reference/cfrank.md)
  (Frank copula)

## See also

[`ddcopula()`](https://janolefi.github.io/RTMBdist/reference/ddcopula.md),
[`dmvcopula()`](https://janolefi.github.io/RTMBdist/reference/dmvcopula.md)

## Examples

``` r
# Normal + Exponential margins with Gaussian copula
x <- c(0.5, 1); y <- c(1, 2)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
dcopula(d1, d2, p1, p2, copula = cgaussian(0.5), log = TRUE)
#> [1] -2.818014 -4.809862

# Normal + Beta margins with Clayton copula
x <- c(0.5, 1); y <- c(0.2, 0.8)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
dcopula(d1, d2, p1, p2, copula = cclayton(2), log = TRUE)
#> [1] -3.8093660 -0.1671136

# Normal + Beta margins with Gumbel copula
x <- c(0.5, 1); y <- c(0.2, 0.4)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
dcopula(d1, d2, p1, p2, copula = cgumbel(1.5), log = TRUE)
#> [1] -1.807264 -1.274899

# Normal + Exponential margins with Frank copula
x <- c(0.5, 1); y <- c(1, 2)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
dcopula(d1, d2, p1, p2, copula = cfrank(2), log = TRUE)
#> [1] -4.585248 -7.325831
```
