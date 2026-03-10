# Multivariate Gaussian copula constructor

Returns a function computing the log density of the multivariate
Gaussian copula, intended to be used with
[`dmvcopula`](https://janolefi.github.io/RTMBdist/reference/dmvcopula.md).

## Usage

``` r
cmvgauss(R)
```

## Arguments

- R:

  Positive definite correlation matrix (unit diagonal)

## Value

Function with matrix argument `U` returning log copula density.

## See also

[`cgmrf()`](https://janolefi.github.io/RTMBdist/reference/cgmrf.md)

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
