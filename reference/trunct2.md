# Truncated t distribution with location and scale

Density, distribution function, quantile function, and random generation
for the truncated t distribution with location `mu` and scale `sigma`.

## Usage

``` r
dtrunct2(x, df, mu = 0, sigma = 1, min = -Inf, max = Inf, log = FALSE)

ptrunct2(q, df, mu = 0, sigma = 1, min = -Inf, max = Inf,
         lower.tail = TRUE, log.p = FALSE)

qtrunct2(p, df, mu = 0, sigma = 1, min = -Inf, max = Inf,
         lower.tail = TRUE, log.p = FALSE)

rtrunct2(n, df, mu = 0, sigma = 1, min = -Inf, max = Inf)
```

## Arguments

- x, q:

  vector of quantiles

- df:

  degrees of freedom parameter, must be positive.

- mu:

  location parameter.

- sigma:

  scale parameter, must be positive.

- min, max:

  truncation bounds.

- log, log.p:

  logical; if `TRUE`, probabilities/densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return.

## Value

`dtrunct2` gives the density, `ptrunct2` gives the distribution
function, `qtrunct2` gives the quantile function, and `rtrunct2`
generates random deviates.

## Details

This implementation of `dtrunct2` allows for automatic differentiation
with `RTMB`.

\$\$f(x;\\\nu,\mu,\sigma,a,b) =
\frac{f_t((x-\mu)/\sigma;\\\nu)/\sigma}{F_t((b-\mu)/\sigma;\\\nu) -
F_t((a-\mu)/\sigma;\\\nu)}, \quad x \in \[a, b\].\$\$

## Examples

``` r
x <- rtrunct2(1, df = 5, mu = 2, sigma = 3, min = -1, max = 5)
d <- dtrunct2(x, df = 5, mu = 2, sigma = 3, min = -1, max = 5)
p <- ptrunct2(x, df = 5, mu = 2, sigma = 3, min = -1, max = 5)
q <- qtrunct2(p, df = 5, mu = 2, sigma = 3, min = -1, max = 5)
```
