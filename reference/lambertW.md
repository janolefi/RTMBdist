# Lambert W function (principal branch)

Solves \\W(x) e^{W(x)} = x\\ for the principal branch \\W_0\\.

## Usage

``` r
lambertW(x)
```

## Arguments

- x:

  vector of evaluation points, \\x \ge -1/e\\.

## Value

The principal branch of the Lambert W function evaluated at `x`.

## Details

This implementation allows for automatic differentiation with `RTMB`.

The value is obtained by Halley iteration, which converges to machine
precision in a handful of steps over the whole domain. For AD, the
function is registered as an atomic operation via
[`ADjoint`](https://rdrr.io/pkg/RTMB/man/ADjoint.html) with the analytic
derivative \$\$W'(x) = \frac{1}{e^{W(x)} (1 + W(x))},\$\$ expressed
through the returned value rather than through \\x\\. Written this way
the derivative is itself an AD-able expression, so derivatives of every
order are available; in particular the third-order derivatives that the
gradient of a Laplace approximation requires.

The principal branch is defined for \\x \ge -1/e\\, with \\W_0(-1/e) =
-1\\ and \\W_0(x) \ge -1\\ throughout. Values below \\-1/e\\ return
`NaN` with a warning. The derivative is infinite at the branch point \\x
= -1/e\\.

## Examples

``` r
lambertW(exp(1)) # 1
#> [1] 1
lambertW(0) # 0
#> [1] 0
x <- c(0.5, 1, 10, 1000)
lambertW(x) * exp(lambertW(x)) - x # ~ 0
#> [1]  0.000000e+00  0.000000e+00 -3.552714e-15 -4.547474e-13
```
