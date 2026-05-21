# Truncated t distribution

Density, distribution function, quantile function, and random generation
for the truncated t distribution.

## Usage

``` r
dtrunct(x, df, min = -Inf, max = Inf, log = FALSE)

ptrunct(q, df, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE)

qtrunct(p, df, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE)

rtrunct(n, df, min = -Inf, max = Inf)
```

## Arguments

- x, q:

  vector of quantiles

- df:

  degrees of freedom parameter, must be positive.

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

`dtrunct` gives the density, `ptrunct` gives the distribution function,
`qtrunct` gives the quantile function, and `rtrunct` generates random
deviates.

## Details

This implementation of `dtrunct` allows for automatic differentiation
with `RTMB`.

\$\$f(x;\\\nu,a,b) = \frac{f_t(x;\\\nu)}{F_t(b;\\\nu) - F_t(a;\\\nu)},
\quad x \in \[a, b\],\$\$ where \\f_t(\cdot;\nu)\\ and
\\F_t(\cdot;\nu)\\ are the Student-\\t\\ PDF and CDF.

## Examples

``` r
x <- rtrunct(1, df = 5, min = -1, max = 5)
d <- dtrunct(x, df = 5, min = -1, max = 5)
p <- ptrunct(x, df = 5, min = -1, max = 5)
q <- qtrunct(p, df = 5, min = -1, max = 5)
```
