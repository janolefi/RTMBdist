# Zero-truncated beta-binomial distribution

Probability mass function, distribution function and random generation
for the zero-truncated beta-binomial distribution.

## Usage

``` r
dztbetabinom(x, size, shape1, shape2, log = FALSE)

pztbetabinom(q, size, shape1, shape2, lower.tail = TRUE, log.p = FALSE)

rztbetabinom(n, size, shape1, shape2)
```

## Arguments

- x:

  integer vector of counts

- size:

  number of trials (zero or more)

- shape1, shape2:

  positive shape parameters of the mixing beta distribution

- log:

  logical; return log-density if TRUE

- q:

  vector of quantiles.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le q\]\\,
  otherwise \\P\[X \> q\]\\.

- log.p:

  logical; if `TRUE`, probabilities are returned on the log scale.

- n:

  number of random values to return.

## Value

`dztbetabinom` gives the probability mass function, `pztbetabinom` gives
the distribution function, and `rztbetabinom` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

By definition, this distribution only has support on the positive
integers (1, ..., n). Any zero-truncated distribution is defined as
\$\$P(X=x \| X\>0) = P(X=x) / (1 - P(X=0)),\$\$ where \\P(X=x)\\ is the
probability mass function of the corresponding untruncated distribution.

The distribution function has no closed form and is computed by summing
the probability mass function over \\0, \ldots, q\\. It is AD-compatible
in the parameters, while `q` and `size` must be numeric data. This is
also what one-step-ahead (OSA) residuals via
`RTMB::`[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
need, so these are supported, e.g. with `method = "cdf"` and
`discrete = TRUE`.

## See also

[betabinom](https://janolefi.github.io/RTMBdist/reference/betabinom.md),
[zibetabinom](https://janolefi.github.io/RTMBdist/reference/zibetabinom.md),
[hbetabinom](https://janolefi.github.io/RTMBdist/reference/hbetabinom.md),
[ztbinom](https://janolefi.github.io/RTMBdist/reference/ztbinom.md)

## Examples

``` r
set.seed(123)
x <- rztbetabinom(5, size = 10, shape1 = 2, shape2 = 3)
d <- dztbetabinom(x, size = 10, shape1 = 2, shape2 = 3)
p <- pztbetabinom(x, size = 10, shape1 = 2, shape2 = 3)
```
