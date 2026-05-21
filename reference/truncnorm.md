# Truncated normal distribution

Density, distribution function, quantile function, and random generation
for the truncated normal distribution.

## Usage

``` r
dtruncnorm(x, mean = 0, sd = 1, min = -Inf, max = Inf, log = FALSE)

ptruncnorm(q, mean = 0, sd = 1, min = -Inf, max = Inf,
           lower.tail = TRUE, log.p = FALSE)

qtruncnorm(p, mean = 0, sd = 1, min = -Inf, max = Inf,
           lower.tail = TRUE, log.p = FALSE)

rtruncnorm(n, mean = 0, sd = 1, min = -Inf, max = Inf)
```

## Arguments

- x, q:

  vector of quantiles

- mean:

  mean parameter, must be positive.

- sd:

  standard deviation parameter, must be positive.

- min, max:

  truncation bounds.

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

`dtruncnorm` gives the density, `ptruncnorm` gives the distribution
function, `qtruncnorm` gives the quantile function, and `rtruncnorm`
generates random deviates.

## Details

This implementation of `dtruncnorm` allows for automatic differentiation
with `RTMB`.

\$\$f(x;\\\mu,\sigma,a,b) =
\frac{\phi((x-\mu)/\sigma)/\sigma}{\Phi((b-\mu)/\sigma) -
\Phi((a-\mu)/\sigma)}, \quad x \in \[a, b\],\$\$ where \\\phi\\ and
\\\Phi\\ are the standard normal PDF and CDF.

## Examples

``` r
x <- rtruncnorm(1, mean = 2, sd = 2, min = -1, max = 5)
d <- dtruncnorm(x, mean = 2, sd = 2, min = -1, max = 5)
p <- ptruncnorm(x, mean = 2, sd = 2, min = -1, max = 5)
q <- qtruncnorm(p, mean = 2, sd = 2, min = -1, max = 5)
```
