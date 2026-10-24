# Generalised beta distribution of the second kind and its special cases

Density, distribution function, quantile function, and random generation
for the generalised beta distribution of the second kind (GB2) and three
of its special cases: the Dagum, Singh-Maddala (Burr XII) and Lomax
(Pareto II) distributions.

## Usage

``` r
dgb2(x, scale = 1, shape1, shape2, shape3, log = FALSE)

pgb2(q, scale = 1, shape1, shape2, shape3, lower.tail = TRUE, log.p = FALSE)

qgb2(p, scale = 1, shape1, shape2, shape3, lower.tail = TRUE, log.p = FALSE)

rgb2(n, scale = 1, shape1, shape2, shape3)

ddagum(x, scale = 1, shape1, shape2, log = FALSE)

pdagum(q, scale = 1, shape1, shape2, lower.tail = TRUE, log.p = FALSE)

qdagum(p, scale = 1, shape1, shape2, lower.tail = TRUE, log.p = FALSE)

rdagum(n, scale = 1, shape1, shape2)

dsinmad(x, scale = 1, shape1, shape3, log = FALSE)

psinmad(q, scale = 1, shape1, shape3, lower.tail = TRUE, log.p = FALSE)

qsinmad(p, scale = 1, shape1, shape3, lower.tail = TRUE, log.p = FALSE)

rsinmad(n, scale = 1, shape1, shape3)

dlomax(x, scale = 1, shape3, log = FALSE)

plomax(q, scale = 1, shape3, lower.tail = TRUE, log.p = FALSE)

qlomax(p, scale = 1, shape3, lower.tail = TRUE, log.p = FALSE)

rlomax(n, scale = 1, shape3)
```

## Arguments

- x, q:

  vector of quantiles.

- scale:

  scale parameter \\b\\, must be positive.

- shape1:

  shape parameter \\a\\, must be positive.

- shape2:

  shape parameter \\p\\, must be positive.

- shape3:

  shape parameter \\q\\, must be positive.

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

`d...` gives the density, `p...` gives the distribution function, `q...`
gives the quantile function, and `r...` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

The GB2 distribution with scale \\b\\ and shape parameters \\a\\, \\p\\
and \\q\\ has density \$\$f(x) = \frac{a\\ x^{ap - 1}}{b^{ap} B(p, q)
\left(1 + (x/b)^a\right)^{p + q}}, \quad x \> 0,\$\$ and its
distribution function is a regularised incomplete beta function,
\$\$F(x) = I_z(p, q), \qquad z = \frac{(x/b)^a}{1 + (x/b)^a}.\$\$ The
special cases fix some of the shapes, and their arguments keep the names
of the shapes that remain:

- Dagum: \\q = 1\\, `ddagum(x, scale, shape1, shape2)`.

- Singh-Maddala, or Burr type XII: \\p = 1\\,
  `dsinmad(x, scale, shape1, shape3)`.

- Lomax, or Pareto type II: \\a = p = 1\\, `dlomax(x, scale, shape3)`.

The log-logistic distribution
[`dllogis`](https://janolefi.github.io/RTMBdist/reference/llogis.md) is
the GB2 with \\p = q = 1\\, and the beta prime distribution
[`dbetaprime`](https://janolefi.github.io/RTMBdist/reference/betaprime.md)
the GB2 with \\a = b = 1\\. The \\k\\-th moment exists for \\-ap \< k \<
aq\\.

The distribution functions are differentiable, so one-step-ahead (OSA)
residuals via
`RTMB::`[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
with `method = "cdf"` are supported.

## References

McDonald, J. B. (1984) Some generalized functions for the size
distribution of income. Econometrica, 52, 647-663.

Kleiber, C. and Kotz, S. (2003) Statistical Size Distributions in
Economics and Actuarial Sciences. Wiley.

## See also

[llogis](https://janolefi.github.io/RTMBdist/reference/llogis.md),
[betaprime](https://janolefi.github.io/RTMBdist/reference/betaprime.md),
[pareto](https://janolefi.github.io/RTMBdist/reference/pareto.md)

## Examples

``` r
set.seed(1)
x <- rgb2(5, scale = 2, shape1 = 3, shape2 = 1.5, shape3 = 0.8)
d <- dgb2(x, 2, 3, 1.5, 0.8)
p <- pgb2(x, 2, 3, 1.5, 0.8)
q <- qgb2(p, 2, 3, 1.5, 0.8)

d <- ddagum(x, 2, 3, 1.5)
d <- dsinmad(x, 2, 3, 0.8)
d <- dlomax(x, 2, 0.8)
```
