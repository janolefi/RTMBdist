# Johnson-Wehrly circular-linear copula constructor

Returns a function computing the log density of the circular-linear
copula of Johnson and Wehrly (1978), intended to be used with
[`dcopula`](https://janolefi.github.io/RTMBdist/reference/dcopula.md) to
join a circular and a linear margin, such as the turning angles and step
lengths of an animal track.

## Usage

``` r
cjw(g, ..., q = 1)
```

## Arguments

- g:

  circular density function with first argument `x` and an argument
  `log`, such as
  [`dvm`](https://janolefi.github.io/RTMBdist/reference/vm.md) or
  [`dwrpcauchy`](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md).

- ...:

  parameters passed to `g`, such as `mu` and `kappa` for
  [`dvm`](https://janolefi.github.io/RTMBdist/reference/vm.md) or `mu`
  and `rho` for
  [`dwrpcauchy`](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md).

- q:

  direction of the dependence, either `1` or `-1`.

## Value

Function of two arguments `(u, v)` returning the log copula density,
with `u` for the circular and `v` for the linear margin.

## Details

The copula density is \$\$c(u, v) = 2\pi \\ g\bigl(2\pi(u - q
v)\bigr),\$\$ where \\u\\ is the distribution function of the circular
margin, \\v\\ that of the linear margin, \\g\\ is a density on the
circle, the *binding density*, and \\q = \pm 1\\. The joint density of
an angle \\\theta\\ and a step length \\s\\ is then \$\$f(\theta, s) =
2\pi \\ g\bigl(2\pi(F_1(\theta) - q F_2(s))\bigr) f_1(\theta)
f_2(s).\$\$ Because \\g\\ is periodic, \\c\\ is a copula for any
circular density and any integer \\q \neq 0\\; only \\q = \pm 1\\ is
allowed here. Any circular density of this package that allows for
automatic differentiation can be used for \\g\\, in particular
[`dvm`](https://janolefi.github.io/RTMBdist/reference/vm.md) and
[`dwrpcauchy`](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md).
Its concentration controls the strength of the dependence, and `q = -1`
reverses its direction.

In
[`dcopula`](https://janolefi.github.io/RTMBdist/reference/dcopula.md),
the circular margin must be the *first* one, i.e. `d1` and `p1` belong
to the angle. Its distribution function can be
[`pwrpcauchy`](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md),
which allows for automatic differentiation. Where the circle is cut open
for this distribution function only shifts \\u\\ by a constant, which
the location of \\g\\ absorbs. The likelihood is therefore the same for
any origin, but the meaning of the location of \\g\\ is not.

The dependence of this copula is a helix: as \\v\\ goes from 0 to 1, the
angle it favours turns once around the circle. The copula is hence not
symmetric in the sign of the angle, \\c(u, v) \neq c(1 - u, v)\\. With a
turning angle margin centred at 0 and the location of \\g\\ at 0, for
example, medium steps are straight, short steps turn one way and long
steps the other. It therefore cannot capture the pattern most common in
movement data, where long steps are straight and short steps turn in
either direction. For this, see
[`cfold`](https://janolefi.github.io/RTMBdist/reference/cfold.md).

Random pairs from the copula are obtained by drawing \\v\\ uniformly and
\\z\\ from \\g\\, and setting \\u = (z / (2\pi) + q v) \bmod 1\\; see
the examples.

## References

Johnson, R. A. and Wehrly, T. E. (1978) Some angular-linear
distributions and related regression models. Journal of the American
Statistical Association, 73, 602-606,
doi:10.1080/01621459.1978.10480062.

Hodel, F. H. and Fieberg, J. R. (2022) Circular-linear copulae for
animal movement data. Methods in Ecology and Evolution, 13,
doi:10.1111/2041-210X.13821.

## See also

[`dcopula()`](https://janolefi.github.io/RTMBdist/reference/dcopula.md),
[`cfold()`](https://janolefi.github.io/RTMBdist/reference/cfold.md),
[vm](https://janolefi.github.io/RTMBdist/reference/vm.md),
[wrpcauchy](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md)

## Examples

``` r
# turning angles with wrapped Cauchy margin, step lengths with Weibull margin
angle <- c(-2, -0.3, 0.1, 1.5); step <- c(0.2, 1.1, 2.4, 0.6)
d1 <- dwrpcauchy(angle, 0, 0.5, log = TRUE); p1 <- pwrpcauchy(angle, 0, 0.5)
d2 <- dweibull(step, 2, 1, log = TRUE); p2 <- pweibull(step, 2, 1)

# von Mises binding density
dcopula(d1, d2, p1, p2, copula = cjw(dvm, mu = 0, kappa = 2), log = TRUE)
#> [1] -2.447020 -3.192045 -7.674264 -4.992124

# wrapped Cauchy binding density, dependence in the other direction
dcopula(d1, d2, p1, p2, copula = cjw(dwrpcauchy, mu = 0, rho = 0.7, q = -1), log = TRUE)
#> [1] -3.326586 -0.436645 -6.666231 -3.140293

# simulation from the joint distribution with von Mises binding density
n <- 1000
v <- runif(n)
z <- rvm(n, mu = 0, kappa = 2)
u <- (z / (2 * pi) + v) %% 1
angle <- qwrpcauchy(u, 0, 0.5)
step <- qweibull(v, 2, 1)
```
