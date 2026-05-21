# Reparameterised gamma distribution

Density, distribution function, quantile function, and random generation
for the gamma distribution reparameterised in terms of mean and standard
deviation.

## Usage

``` r
dgamma2(x, mean = 1, sd = 1, log = FALSE)

pgamma2(q, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE)

qgamma2(p, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE)

rgamma2(n, mean = 1, sd = 1)
```

## Arguments

- x, q:

  vector of quantiles

- mean:

  mean parameter, must be positive.

- sd:

  standard deviation parameter, must be positive.

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

`dgamma2` gives the density, `pgamma2` gives the distribution function,
`qgamma2` gives the quantile function, and `rgamma2` generates random
deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

Reparameterises the gamma distribution via \\\text{shape} = \mu^2/s^2\\
and \\\text{scale} = s^2/\mu\\: \$\$f(x;\\\mu, s) = \frac{x^{\mu^2/s^2 -
1} \exp(-x s^2/\mu^2)}{(s^2/\mu)^{\mu^2/s^2}\\\Gamma(\mu^2/s^2)}, \quad
x \> 0.\$\$

## Examples

``` r
x <- rgamma2(1)
d <- dgamma2(x)
p <- pgamma2(x)
q <- qgamma2(p)
```
