# Log-logistic distribution

Density, distribution function, quantile function, and random generation
for the log-logistic distribution.

## Usage

``` r
dllogis(x, alpha = 1, beta = 1, log = FALSE)

pllogis(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)

qllogis(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)

rllogis(n, alpha = 1, beta = 1)
```

## Arguments

- x, q:

  vector of quantiles (\\x \> 0\\).

- alpha:

  scale parameter (\\\alpha \> 0\\); equal to the median.

- beta:

  shape parameter (\\\beta \> 0\\); controls tail heaviness.

- log:

  logical; if `TRUE`, densities are returned on the log scale.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- log.p:

  logical; if `TRUE`, probabilities are returned on the log scale.

- p:

  vector of probabilities.

- n:

  number of random values to return.

## Value

`dllogis` gives the density, `pllogis` gives the distribution function,
`qllogis` gives the quantile function, and `rllogis` generates random
deviates.

## Details

The log-logistic distribution has density \$\$f(x;\\\alpha,\beta) =
\frac{(\beta/\alpha)\\(x/\alpha)^{\beta-1}}{\bigl(1 +
(x/\alpha)^\beta\bigr)^2}, \quad x \> 0,\$\$ where \\\alpha \> 0\\ is
the scale parameter and \\\beta \> 0\\ is the shape parameter. The scale
parameter equals the median. Larger \\\beta\\ gives lighter tails; the
distribution has a finite mean only when \\\beta \> 1\\ and a finite
variance only when \\\beta \> 2\\.

The log-logistic arises naturally as the distribution of \\X = e^Y\\
where \\Y \sim \mathrm{Logistic}(\log\alpha,\\ 1/\beta)\\, which yields
the numerically convenient log-density \$\$\log f(x) = \log\beta - \log
x + \log f\_{\mathrm{logistic}}\\\bigl(\beta\log(x/\alpha)\bigr).\$\$

The CDF is \$\$F(x;\\\alpha,\beta) = \frac{1}{1 +
(x/\alpha)^{-\beta}},\$\$ and the quantile function is
\$\$Q(p;\\\alpha,\beta) = \alpha
\left(\frac{p}{1-p}\right)^{1/\beta}.\$\$

`dllogis` allows for automatic differentiation with `RTMB`.

## Examples

``` r
x <- rllogis(5, alpha = 1, beta = 2)
dllogis(x, alpha = 1, beta = 2)
#> [1] 0.06015175 0.64946566 0.38456522 0.63102299 0.01896198
pllogis(c(0.5, 1, 2), alpha = 1, beta = 2)
#> [1] 0.2 0.5 0.8
qllogis(c(0.25, 0.5, 0.75), alpha = 1, beta = 2)
#> [1] 0.5773503 1.0000000 1.7320508
```
