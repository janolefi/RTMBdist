# Half-t distribution

Density, distribution function, quantile function, and random generation
for the half-t distribution.

## Usage

``` r
dhalft(x, df, sigma = 1, log = FALSE)

phalft(q, df, sigma = 1, lower.tail = TRUE, log.p = FALSE)

qhalft(p, df, sigma = 1, lower.tail = TRUE, log.p = FALSE)

rhalft(n, df, sigma = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- df:

  degrees of freedom, must be positive.

- sigma:

  scale parameter, must be positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities.

- n:

  number of random values to return.

## Value

`dhalft` gives the density, `phalft` gives the distribution function,
`qhalft` gives the quantile function, and `rhalft` generates random
deviates.

## Details

`dhalft` and `phalft` allow for automatic differentiation with `RTMB`,
with respect to `df` as well as `sigma`.

The half-t is the distribution of \\\|Y\|\\ for \\Y = \sigma T\\ and
\\T\\ Student t on \\\nu\\ degrees of freedom: \$\$f(x;\\\nu,\sigma) =
\frac{2}{\sigma} f_T(x/\sigma;\\ \nu), \quad x \ge 0,\$\$ with
distribution function \\F(x) = 2 F_T(x/\sigma;\\ \nu) - 1\\.

Together with the
[half-Cauchy](https://janolefi.github.io/RTMBdist/reference/halfcauchy.md),
which is the case `df = 1`, this is the standard weakly informative
prior for the standard deviation of a hierarchical model, recommended by
Gelman (2006) in place of the inverse-gamma. The degrees of freedom set
how heavy the tail is: small `df` leaves large values essentially
unpenalised, and as `df` grows the distribution approaches the
half-normal, which is the [folded
normal](https://janolefi.github.io/RTMBdist/reference/foldnorm.md) with
`mu = 0`.

The mean \$\$E(X) =
\frac{2\sigma\sqrt{\nu}\\\Gamma\bigl(\tfrac{\nu+1}{2}\bigr)}{\sqrt{\pi}\\(\nu -
1)\\\Gamma\bigl(\tfrac{\nu}{2}\bigr)}\$\$ exists only for \\\nu \> 1\\
and the variance \\\sigma^2 \nu / (\nu - 2) - E(X)^2\\ only for \\\nu \>
2\\.

## References

Gelman, A. (2006) Prior distributions for variance parameters in
hierarchical models. Bayesian Analysis, 1, 515-534,
doi:10.1214/06-BA117A.

## See also

[halfcauchy](https://janolefi.github.io/RTMBdist/reference/halfcauchy.md),
[foldnorm](https://janolefi.github.io/RTMBdist/reference/foldnorm.md),
[trunct](https://janolefi.github.io/RTMBdist/reference/trunct.md),
[t2](https://janolefi.github.io/RTMBdist/reference/t2.md)

## Examples

``` r
set.seed(123)
x <- rhalft(5, df = 3, sigma = 2)
d <- dhalft(x, df = 3, sigma = 2)
p <- phalft(x, df = 3, sigma = 2)
q <- qhalft(p, df = 3, sigma = 2)
```
