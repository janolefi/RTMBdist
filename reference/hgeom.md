# Hurdle geometric distribution

Probability mass function, distribution function, and random generation
for the hurdle (zero-altered) geometric distribution.

## Usage

``` r
dhgeom(x, prob, zeroprob = 0.5, log = FALSE)

phgeom(q, prob, zeroprob = 0.5, lower.tail = TRUE, log.p = FALSE)

rhgeom(n, prob, zeroprob = 0.5)
```

## Arguments

- x, q:

  integer vector of counts

- prob:

  probability of success in each trial, in (0,1)

- zeroprob:

  probability of a zero, between 0 and 1

- log, log.p:

  logical; return log-density if TRUE

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dhgeom` gives the probability mass function, `phgeom` gives the
distribution function, and `rhgeom` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

A hurdle distribution models the zeros and the positive counts as two
separate processes: the probability of a zero is a free parameter, and
the positive counts follow the corresponding zero-truncated
distribution. Writing \\p_0\\ for `zeroprob`, \$\$P(X = 0) = p_0, \qquad
P(X = x) = (1 - p_0)\\\frac{P\_{\mathrm{Geom}}(x;\\\pi)}{1 - \pi_0},
\quad x = 1, 2, \ldots\$\$ where \\\pi_0 = P\_{\mathrm{Geom}}(0;\\\pi) =
\pi\\ is the probability of a zero under the ordinary geometric.

Unlike zero-inflation, which can only add zeros to those the geometric
already produces, `zeroprob` here is exactly the probability of a zero
and may be larger *or* smaller than \\\pi_0\\. The two coincide with the
ordinary geometric when `zeroprob` equals \\\pi_0\\.

## References

Mullahy, J. (1986) Specification and testing of some modified count data
models. Journal of Econometrics, 33, 341-365.

## See also

[zigeom](https://janolefi.github.io/RTMBdist/reference/zigeom.md),
[ztgeom](https://janolefi.github.io/RTMBdist/reference/ztgeom.md),
[hnbinom](https://janolefi.github.io/RTMBdist/reference/hnbinom.md),
[hpois](https://janolefi.github.io/RTMBdist/reference/hpois.md)

## Examples

``` r
set.seed(123)
x <- rhgeom(5, prob = 0.3, zeroprob = 0.4)
d <- dhgeom(x, prob = 0.3, zeroprob = 0.4)
p <- phgeom(x, prob = 0.3, zeroprob = 0.4)
```
