# von Mises-Fisher distribution

Density, distribution function, and random generation for the von
Mises-Fisher distribution.

## Usage

``` r
dvmf(x, mu, kappa, log = FALSE)

rvmf(n, mu, kappa)
```

## Arguments

- x:

  unit vector or matrix (with each row being a unit vector) of
  evaluation points

- mu:

  unit mean vector

- kappa:

  non-negative numeric value for the concentration parameter of the
  distribution.

- log:

  logical; if `TRUE`, densities are returned on the log scale.

- n:

  number of random values to return.

## Value

`dvmf` gives the density and `rvm` generates random deviates.

## Details

This implementation of `dvmf` allows for automatic differentiation with
`RTMB`. `rvmf` is a reparameterised import from
[`movMF::rmovMF`](https://rdrr.io/pkg/movMF/man/movMF_distribution.html).

\$\$f(\mathbf{x};\\\boldsymbol{\mu},\kappa) =
C_p(\kappa)\\\exp(\kappa\\\boldsymbol{\mu}^\top\mathbf{x}),\$\$ where
\\C_p(\kappa) = \kappa^{p/2-1} / ((2\pi)^{p/2} I\_{p/2-1}(\kappa))\\ is
the normalising constant and \\I\_\nu\\ is the modified Bessel function
of the first kind of order \\\nu\\.

## Examples

``` r
set.seed(123)
# single parameter set
mu <- rep(1, 3) / sqrt(3)
kappa <- 4
x <- rvmf(1, mu, kappa)
d <- dvmf(x, mu, kappa)

# vectorised over parameters
mu <- matrix(mu, nrow = 1)
mu <- mu[rep(1,10), ]
kappa <- rep(kappa, 10)
x <- rvmf(10, mu, kappa)
d <- dvmf(x, mu, kappa)
```
