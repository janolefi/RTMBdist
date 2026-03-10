# Joint probability under a discrete bivariate copula

Computes the joint probability mass function of two **discrete** margins
combined with a copula CDF.

## Usage

``` r
ddcopula(d1, d2, p1, p2, Copula, log = FALSE)
```

## Arguments

- d1, d2:

  Marginal p.m.f. values at the observed points, **not** on log-scale.

- p1, p2:

  Marginal CDF values at the observed points.

- Copula:

  A function of two arguments returning the copula CDF.

- log:

  Logical; if `TRUE`, return the log joint density. In this case, `d1`
  and `d2` must be on the log scale.

## Value

Joint probability (or log-probability) under chosen copula

## Details

The joint probability mass function for two discrete margins is \$\$
\Pr(Y_1 = y_1, Y_2 = y_2) = C(F_1(y_1), F_2(y_2)) - C(F_1(y_1-1),
F_2(y_2)) - C(F_1(y_1), F_2(y_2-1)) + C(F_1(y_1-1), F_2(y_2-1)), \$\$
where \\F_i\\ are the marginal CDFs, and \\C\\ is the copula CDF.

Available copula CDF constructors are:

- [`Cclayton`](https://janolefi.github.io/RTMBdist/reference/cclayton.md)
  (Clayton copula)

- [`Cgumbel`](https://janolefi.github.io/RTMBdist/reference/cgumbel.md)
  (Gumbel copula)

- [`Cfrank`](https://janolefi.github.io/RTMBdist/reference/cfrank.md)
  (Frank copula)

## See also

[`dcopula()`](https://janolefi.github.io/RTMBdist/reference/dcopula.md),
[`dmvcopula()`](https://janolefi.github.io/RTMBdist/reference/dmvcopula.md)

## Examples

``` r
x <- c(3,5); y <- c(2,4)
d1 <- dpois(x, 4); d2 <- dpois(y, 3)
p1 <- ppois(x, 4); p2 <- ppois(y, 3)
ddcopula(d1, d2, p1, p2, Copula = Cclayton(2), log = FALSE)
#> [1] 0.07616878 0.04042844
```
