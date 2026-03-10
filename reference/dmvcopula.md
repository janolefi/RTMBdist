# Joint density under a multivariate copula

Computes the joint density (or log-density) of a distribution
constructed from any number of arbitrary margins combined with a
specified copula.

## Usage

``` r
dmvcopula(D, P, copula = cmvgauss(diag(ncol(D))), log = FALSE)
```

## Arguments

- D:

  Matrix of marginal density values of with rows corresponding to
  observations and columns corresponding to dimensions. If `log = TRUE`,
  supply the log-densities. If `log = FALSE`, supply the raw densities

- P:

  Matrix of marginal CDF values of the same dimension as `D`. Need not
  be supplied on log scale.

- copula:

  A function of a matrix argument `U` returning the log copula density
  \\\log c(u_1, ... u_d)\\. The columns of `U` correspond to dimensions.
  You can either construct this yourself or use the copula constructors
  available (see details)

- log:

  Logical; if `TRUE`, return the log joint density. In this case, `D`
  must be on the log scale.

## Value

Joint density (or log-density) under the chosen copula.

## Details

The joint density is \$\$f(x_1, \dots, x_d) = c(F_1(x_1), \dots,
F_d(x_d)) \\ f_1(x_1) \dots f_d(x_d),\$\$ where \\F_i\\ are the marginal
CDFs, \\f_i\\ are the marginal densities, and \\c\\ is the copula
density.

The marginal densities `d_1, ..., d_d` and CDFs `p_1, ..., p_d` must be
differentiable for automatic differentiation (AD) to work.

Available multivariate copula constructors are:

- [`cmvgauss`](https://janolefi.github.io/RTMBdist/reference/cmvgauss.md)
  (Multivariate Gaussian copula)

- [`cgmrf`](https://janolefi.github.io/RTMBdist/reference/cgmrf.md)
  (Multivariate Gaussian copula parameterised by precision (inverse
  correlation) matrix)

## See also

[`dcopula()`](https://janolefi.github.io/RTMBdist/reference/dcopula.md),
[`ddcopula()`](https://janolefi.github.io/RTMBdist/reference/ddcopula.md)

## Examples

``` r
x <- c(0.5, 1); y <- c(1, 2); z <- c(0.2, 0.8)
d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE); d3 <- dbeta(z, 2, 1, log = TRUE)
p1 <- pnorm(x, 1); p2 <- pexp(y, 2); p3 <- pbeta(z, 2, 1)
R <- matrix(c(1,0.5,0.3,0.5,1,0.4,0.3,0.4,1), nrow = 3)

### Multivariate Gaussian copula
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
