# Half-Cauchy distribution

Density, distribution function, quantile function, and random generation
for the half-Cauchy distribution.

## Usage

``` r
dhalfcauchy(x, sigma = 1, log = FALSE)

phalfcauchy(q, sigma = 1, lower.tail = TRUE, log.p = FALSE)

qhalfcauchy(p, sigma = 1, lower.tail = TRUE, log.p = FALSE)

rhalfcauchy(n, sigma = 1)
```

## Arguments

- x, q:

  vector of quantiles.

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

`dhalfcauchy` gives the density, `phalfcauchy` gives the distribution
function, `qhalfcauchy` gives the quantile function, and `rhalfcauchy`
generates random deviates.

## Details

`dhalfcauchy` and `phalfcauchy` allow for automatic differentiation with
`RTMB`.

The half-Cauchy is the distribution of \\\|Y\|\\ for \\Y\\ Cauchy with
scale \\\sigma\\: \$\$f(x;\\\sigma) = \frac{2}{\pi\sigma\bigl(1 +
(x/\sigma)^2\bigr)}, \quad x \ge 0,\$\$ with distribution function
\\F(x) = \frac{2}{\pi}\arctan(x/\sigma)\\.

It is the standard weakly informative prior for the standard deviation
of a hierarchical model, recommended by Gelman (2006) in place of the
inverse-gamma: the density is flat and non-zero at the origin, so it
does not force the variance component away from zero, while the tail is
heavy enough to leave large values unpenalised. That makes it a natural
companion to models fitted by the Laplace approximation, where the
variance components are exactly the parameters at issue.

It has no moments of any order. It is the
[half-t](https://janolefi.github.io/RTMBdist/reference/halft.md) with
`df = 1`, and the [folded
normal](https://janolefi.github.io/RTMBdist/reference/foldnorm.md) with
`mu = 0` is the corresponding half-normal, which the half-t approaches
as `df` grows.

## References

Gelman, A. (2006) Prior distributions for variance parameters in
hierarchical models. Bayesian Analysis, 1, 515-534,
doi:10.1214/06-BA117A.

## See also

[halft](https://janolefi.github.io/RTMBdist/reference/halft.md),
[foldnorm](https://janolefi.github.io/RTMBdist/reference/foldnorm.md),
[trunct](https://janolefi.github.io/RTMBdist/reference/trunct.md)

## Examples

``` r
set.seed(123)
x <- rhalfcauchy(5, sigma = 2)
d <- dhalfcauchy(x, sigma = 2)
p <- phalfcauchy(x, sigma = 2)
q <- qhalfcauchy(p, sigma = 2)
```
