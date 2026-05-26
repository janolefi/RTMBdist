# Gompertz distribution

Density, distribution function, quantile function, and random generation
for the Gompertz distribution.

## Usage

``` r
dgompertz(x, eta = 1, b = 1, log = FALSE)

pgompertz(q, eta = 1, b = 1, lower.tail = TRUE, log.p = FALSE)

qgompertz(p, eta = 1, b = 1, lower.tail = TRUE, log.p = FALSE)

rgompertz(n, eta = 1, b = 1)
```

## Arguments

- x, q:

  vector of quantiles (non-negative)

- eta:

  shape parameter, must be positive

- b:

  rate parameter, must be positive

- log, log.p:

  logical; if `TRUE`, probabilities/densities are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dgompertz` gives the density, `pgompertz` gives the distribution
function, `qgompertz` gives the quantile function, and `rgompertz`
generates random deviates.

## Details

The Gompertz distribution with shape \\\eta \> 0\\ and rate \\b \> 0\\
has density \$\$f(x;\\\eta,b) =
b\eta\\e^{bx}\exp\\\bigl(-\eta(e^{bx}-1)\bigr), \quad x \geq 0,\$\$ with
CDF \\F(x) = 1 - \exp(-\eta(e^{bx}-1))\\ and quantile function \\Q(p) =
\log(1 - \log(1-p)/\eta)\\/\\b\\.

## References

<https://en.wikipedia.org/wiki/Gompertz_distribution>

## Examples

``` r
x <- rgompertz(1, eta = 1, b = 1)
d <- dgompertz(x, eta = 1, b = 1)
p <- pgompertz(x, eta = 1, b = 1)
q <- qgompertz(p, eta = 1, b = 1)
```
