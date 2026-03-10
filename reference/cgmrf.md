# Multivariate Gaussian copula constructor parameterised by inverse correlation matrix

Returns a function computing the log density of the multivariate
Gaussian copula, parameterised by the inverse correlation matrix.

## Usage

``` r
cgmrf(Q)
```

## Arguments

- Q:

  Inverse of a positive definite correlation matrix with unit diagonal.
  Can either be sparse or dense matrix.

## Value

Function with matrix argument `U` returning log copula density.

## Details

**Caution:** Parameterising the inverse correlation directly is
difficult, as inverting it needs to yield a positive definite matrix
with **unit diagonal**. Hence we still advise parameterising the
correaltion matrix `R` and computing its inverse. This function is
useful when you need access to the precision (i.e. inverse correlation)
in your likelihood function.

## See also

[`cmvgauss()`](https://janolefi.github.io/RTMBdist/reference/cmvgauss.md)

## Examples

``` r
x <- c(0.5, 1); y <- c(1, 2); z <- c(0.2, 0.8)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE); d3 <- dbeta(z, 2, 1, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pexp(y, 2); p3 <- pbeta(z, 2, 1)
R <- matrix(c(1,0.5,0.3,0.5,1,0.4,0.3,0.4,1), nrow = 3)

## Based on correlation matrix
dmvcopula(cbind(d1, d2, d3), cbind(p1, p2, p3), copula = cmvgauss(R), log = TRUE)
#> [1] -4.651470 -4.249599

## Based on precision matrix
Q <- solve(R)
dmvcopula(cbind(d1, d2, d3), cbind(p1, p2, p3), copula = cgmrf(Q), log = TRUE)
#> [1] -4.651470 -4.249599

## Parameterisation inside a model
# using RTMB::unstructured to get a valid correlation matrix
library(RTMB)
d <- 5 # dimension
cor_func <- unstructured(d)
npar <- length(cor_func$parms())
R <- cor_func$corr(rep(0.1, npar))
```
