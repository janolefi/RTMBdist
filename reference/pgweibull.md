# Power generalized Weibull distribution

Survival, hazard, cumulative distribution, density, quantile and
sampling function for the power generalized Weibull (PgW) distribution
with parameters `scale`, `shape` and `powershape`.

## Usage

``` r
spgweibull(q, scale = 1, shape = 1, powershape = 1, log = FALSE)

hpgweibull(x, scale = 1, shape = 1, powershape = 1, log = FALSE)

ppgweibull(q, scale = 1, shape = 1, powershape = 1,
           lower.tail = TRUE, log.p = FALSE)

dpgweibull(x, scale = 1, shape = 1, powershape = 1, log = FALSE)

qpgweibull(p, scale = 1, shape = 1, powershape = 1)

rpgweibull(n, scale = 1, shape = 1, powershape = 1)
```

## Arguments

- scale:

  positive scale parameter

- shape:

  positive shape parameter

- powershape:

  positive power shape parameter

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- x, q:

  vector of quantiles

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of observations

## Value

`dpgweibull` gives the density, `ppgweibull` gives the distribution
function, `qpgweibull` gives the quantile function, and `rpgweibull`
generates random deviates. `spgweibull` gives the survival function and
`hpgweibull` gives the hazard function.

## Details

The survival function of the PgW distribution is: \$\$ S(x) = \exp
\left\\ 1 - \left\[ 1 +
\left(\frac{x}{\theta}\right)^{\nu}\right\]^{\frac{1}{\gamma}} \right\\.
\$\$ The hazard function is \$\$ \frac{\nu}{\gamma\theta^{\nu}}\cdot
x^{\nu-1}\cdot \left\[ 1 +
\left(\frac{x}{\theta}\right)^{\nu}\right\]^{\frac{1}{\gamma-1}} \$\$
The cumulative distribution function is then \\F(x) = 1 - S(x)\\ and the
density function is \\S(x)\cdot h(x)\\.

If both shape parameters equal 1, the PgW distribution reduces to the
exponential distribution (see
[`dexp`](https://rdrr.io/r/stats/Exponential.html)) with \\\texttt{rate}
= 1/\texttt{scale}\\ If the power shape parameter equals 1, the PgW
distribution simplifies to the Weibull distribution (see
[`dweibull`](https://rdrr.io/r/stats/Weibull.html)) with the same
parametrization.

## Examples

``` r
x <- rpgweibull(1, 2, 2, 3)
d <- dpgweibull(x, 2, 2, 3)
p <- ppgweibull(x, 2, 2, 3)
q <- qpgweibull(p, 2, 2, 3)
s <- spgweibull(x, 2, 2, 3)
h <- hpgweibull(x, 2, 2, 3)
```
