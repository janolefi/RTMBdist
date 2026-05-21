# Generalised Poisson distribution

Probability mass function, distribution function, and random generation
for the generalised Poisson distribution.

## Usage

``` r
dgenpois(x, lambda = 1, phi = 1, log = FALSE)

pgenpois(q, lambda = 1, phi = 1, lower.tail = TRUE, log.p = FALSE)

qgenpois(p, lambda = 1, phi = 1,
         lower.tail = TRUE, log.p = FALSE, max.value = 10000)

rgenpois(n, lambda = 1, phi = 1, max.value = 10000)
```

## Arguments

- x, q:

  integer vector of counts

- lambda:

  vector of positive means

- phi:

  vector of non-negative dispersion parameters

- log, log.p:

  logical; return log-density if TRUE

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- max.value:

  a constant, set to the default value of 10000 for how far the
  algorithm should look for `q`.

- n:

  number of random values to return.

## Value

`dgenpois` gives the probability mass function, `pgenpois` gives the
distribution function, `qgenpois` gives the quantile function, and
`rgenpois` generates random deviates.

## Details

This implementation of `dgenpois` allows for automatic differentiation
with `RTMB`. The other functions are imported from
[`gamlss.dist::GPO`](https://rdrr.io/pkg/gamlss.dist/man/GPO.html).

The distribution has mean \\\lambda\\ and variance \\\lambda(1 + \phi
\lambda)^2\\. For \\\phi = 0\\ it reduces to the Poisson distribution,
however \\\phi\\ must be strictly positive here.

\$\$P(X = x;\\\lambda,\phi) = \frac{\lambda\\(1+\phi
x)^{x-1}\\e^{-\lambda(1+\phi x)/(1+\phi\lambda)}}{(1+\phi\lambda)^x\\
x!}, \quad x = 0, 1, 2, \ldots\$\$

## Examples

``` r
set.seed(123)
x <- rgenpois(1, 2, 3)
d <- dgenpois(x, 2, 3)
p <- pgenpois(x, 2, 3)
q <- qgenpois(p, 2, 3)
```
