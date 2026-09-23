# Folded circular-linear copula constructor

Turns a copula for two linear variables into a circular-linear copula
that is symmetric in the sign of the angle, intended to be used with
[`dcopula`](https://janolefi.github.io/RTMBdist/reference/dcopula.md) to
join a turning angle and a step length.

## Usage

``` r
cfold(copula)
```

## Arguments

- copula:

  function of two arguments `(u, v)` returning a log copula density for
  two linear variables, e.g. `cgaussian(0.5)`.

## Value

Function of two arguments `(u, v)` returning the log copula density,
with `u` for the circular and `v` for the linear margin.

## Details

The copula density is \$\$c(u, v) = c_0\bigl(1 - \|2u - 1\|, \\
v\bigr),\$\$ where \\u\\ is the distribution function of the circular
margin, \\v\\ that of the linear margin and \\c_0\\ is the density of
the copula passed as `copula`, e.g.
[`cgaussian`](https://janolefi.github.io/RTMBdist/reference/cgaussian.md),
[`cclayton`](https://janolefi.github.io/RTMBdist/reference/cclayton.md),
[`cgumbel`](https://janolefi.github.io/RTMBdist/reference/cgumbel.md) or
[`cfrank`](https://janolefi.github.io/RTMBdist/reference/cfrank.md).
This is a copula for any \\c_0\\. It is the rectangular patchwork copula
of Hodel and Fieberg (2022), with the linear copula in the rectangle \\u
\le 1/2\\ and its mirror image in the rectangle \\u \> 1/2\\.

The map \\u \mapsto 1 - \|2u - 1\|\\ folds the circle at \\u = 1/2\\.
For a circular margin that is symmetric about its mean direction and
whose distribution function is cut open at the antipode of that mean
direction, as by default in
[`pwrpcauchy`](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md)
and [`pvm`](https://janolefi.github.io/RTMBdist/reference/vm.md), \\u =
1/2\\ is the mean direction. \\1 - \|2u - 1\|\\ is then the distribution
function of the angular distance to the mean direction, counted from the
antipode: it is 1 for angles at the mean direction and 0 for angles
opposite to it. In words, \\c_0\\ links the straightness of a step to
its length. A copula with positive dependence, such as `cgaussian(rho)`
with `rho > 0`, makes long steps straight and lets short steps turn in
either direction, the pattern most common in movement data. Unlike
[`cjw`](https://janolefi.github.io/RTMBdist/reference/cjw.md), the
resulting copula is symmetric, \\c(u, v) = c(1 - u, v)\\.

In
[`dcopula`](https://janolefi.github.io/RTMBdist/reference/dcopula.md),
the circular margin must be the *first* one, i.e. `d1` and `p1` belong
to the angle.

## References

Hodel, F. H. and Fieberg, J. R. (2022) Circular-linear copulae for
animal movement data. Methods in Ecology and Evolution, 13,
doi:10.1111/2041-210X.13821.

Durante, F., Saminger-Platz, S. and Sarkoci, P. (2009) Rectangular
patchwork for bivariate copulas and tail dependence. Communications in
Statistics - Theory and Methods, 38, 2515-2527,
doi:10.1080/03610920802571203.

## See also

[`dcopula()`](https://janolefi.github.io/RTMBdist/reference/dcopula.md),
[`cjw()`](https://janolefi.github.io/RTMBdist/reference/cjw.md),
[vm](https://janolefi.github.io/RTMBdist/reference/vm.md),
[wrpcauchy](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md)

## Examples

``` r
# turning angles with von Mises margin, step lengths with Weibull margin
angle <- c(-2, -0.3, 0.1, 1.5); step <- c(0.2, 1.1, 2.4, 0.6)
d1 <- dvm(angle, 0, 2, log = TRUE); p1 <- pvm(angle, 0, 2)
d2 <- dweibull(step, 2, 1, log = TRUE); p2 <- pweibull(step, 2, 1)
dcopula(d1, d2, p1, p2, copula = cfold(cgaussian(0.5)), log = TRUE)
#> [1] -3.2413130 -0.9366028 -3.9252773 -2.4386354

# wrapped Cauchy margin, which also allows for automatic differentiation
d1 <- dwrpcauchy(angle, 0, 0.5, log = TRUE); p1 <- pwrpcauchy(angle, 0, 0.5)
dcopula(d1, d2, p1, p2, copula = cfold(cclayton(2)), log = TRUE)
#> [1] -3.1479573 -0.8881331 -4.0546310 -1.7932076

# the copula is symmetric in the sign of the angle
cop <- cfold(cclayton(2))
cop(c(0.2, 0.8), 0.9)
#> [1] -0.5099969 -0.5099969
```
