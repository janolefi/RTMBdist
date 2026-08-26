# Reparameterised Bell distribution

Probability mass function, distribution function, quantile function, and
random generation for the Bell distribution reparameterised in terms of
its mean.

## Usage

``` r
dbell2(x, mu, log = FALSE)

pbell2(q, mu, lower.tail = TRUE, log.p = FALSE)

qbell2(p, mu, lower.tail = TRUE, log.p = FALSE)

rbell2(n, mu)
```

## Arguments

- x, q:

  integer vector of counts

- mu:

  vector of positive means

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return.

## Value

`dbell2` gives the probability mass function, `pbell2` gives the
distribution function, `qbell2` gives the quantile function, and
`rbell2` generates random deviates.

## Details

This implementation of `dbell2` and `pbell2` allows for automatic
differentiation with `RTMB` with respect to `mu`.

The Bell distribution has mean \\\mu = \theta e^{\theta}\\, which is a
bijection from \\\theta \> 0\\ to \\\mu \> 0\\ and is inverted by the
principal branch of the Lambert W function, \$\$\theta = W(\mu).\$\$
Every positive mean therefore corresponds to exactly one \\\theta\\. All
four functions simply apply this transformation and hand over to their
[`bell`](https://janolefi.github.io/RTMBdist/reference/bell.md)
counterparts.

In this parameterisation the variance is \\\mu (1 + W(\mu))\\, so the
distribution is always overdispersed relative to the Poisson
distribution, but the degree of overdispersion is determined by the mean
rather than by a free parameter.

[`lambertW`](https://janolefi.github.io/RTMBdist/reference/lambertW.md)
is AD-compatible to arbitrary order, so `mu` may be a parameter of a
model fitted by Laplace approximation.

## References

Castellares, F., Ferrari, S. L. P., and Lemonte, A. J. (2018). On the
Bell distribution and its associated regression model for count data.
*Applied Mathematical Modelling* 56, 172-185.
doi:10.1016/j.apm.2017.12.014

## See also

[`bell`](https://janolefi.github.io/RTMBdist/reference/bell.md) for the
natural parameterisation.

## Examples

``` r
set.seed(123)
x <- rbell2(1, 3)
d <- dbell2(x, 3)
p <- pbell2(x, 3)
q <- qbell2(p, 3)

# the two parameterisations agree
all.equal(dbell2(0:5, 3), dbell(0:5, lambertW(3)))
#> [1] TRUE
```
