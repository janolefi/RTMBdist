# wrapped Cauchy distribution

Density and random generation for the wrapped Cauchy distribution.

## Usage

``` r
dwrpcauchy(x, mu = 0, rho, log = FALSE)

rwrpcauchy(n, mu = 0, rho, wrap = TRUE)
```

## Arguments

- x:

  vector of angles measured in radians at which to evaluate the density
  function.

- mu:

  mean direction of the distribution measured in radians.

- rho:

  concentration parameter of the distribution, must be in the interval
  from 0 to 1.

- log:

  logical; if `TRUE`, densities are returned on the log scale.

- n:

  number of random values to return.

- wrap:

  logical; if `TRUE`, generated angles are wrapped to the interval from
  -pi to pi.

## Value

`dwrpcauchy` gives the density and `rwrpcauchy` generates random
deviates.

## Details

This implementation of `dwrpcauchy` allows for automatic differentiation
with `RTMB`. `rwrpcauchy` is simply a wrapper for
`rwrappedcauchy`imported from `circular`.

\$\$f(x;\\\mu,\rho) = \frac{1}{2\pi}\cdot\frac{1-\rho^2}{1 + \rho^2 -
2\rho\cos(x-\mu)}.\$\$

## Examples

``` r
set.seed(1)
x <- rwrpcauchy(10, 0, 0.5)
d <- dwrpcauchy(x, 0, 0.5)
```
