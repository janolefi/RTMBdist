# AD-compatible geometric distribution

Density and distribution function for the geometric distribution,
written so that they can be taped by `RTMB`.

## Usage

``` r
dgeom.ad(x, prob, log = FALSE)

pgeom.ad(q, prob, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- x, q:

  integer vector of counts

- prob:

  probability of success in each trial, in (0,1\]

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

## Value

`dgeom.ad` gives the probability mass function and `pgeom.ad` gives the
distribution function.

## Details

`stats` already provides the geometric distribution, but its versions
cannot be differentiated. These are AD-compatible replacements, reached
automatically whenever an argument is an AD variable, so
[`stats::dgeom`](https://rdrr.io/r/stats/Geometric.html) and
[`stats::pgeom`](https://rdrr.io/r/stats/Geometric.html) are left
untouched for ordinary use.

The parameterisation is the same as in `stats`: \\X\\ is the number of
failures before the first success, so \$\$P(X = x;\\\pi) = \pi\\(1 -
\pi)^{x}, \quad x = 0, 1, 2, \ldots\$\$ The density is obtained from the
negative binomial with `size = 1`, for which `RTMB` provides an AD
method; the distribution function \\1 - (1-\pi)^{x+1}\\ is elementary.

## See also

[zigeom](https://janolefi.github.io/RTMBdist/reference/zigeom.md),
[ztgeom](https://janolefi.github.io/RTMBdist/reference/ztgeom.md),
[hgeom](https://janolefi.github.io/RTMBdist/reference/hgeom.md)

## Examples

``` r
dgeom.ad(0:5, prob = 0.3)
#> [1] 0.300000 0.210000 0.147000 0.102900 0.072030 0.050421
pgeom.ad(0:5, prob = 0.3)
#> [1] 0.300000 0.510000 0.657000 0.759900 0.831930 0.882351
```
