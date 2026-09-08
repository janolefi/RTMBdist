# Zero-truncated beta-binomial distribution

Probability mass function and random generation for the zero-truncated
beta-binomial distribution.

## Usage

``` r
dztbetabinom(x, size, shape1, shape2, log = FALSE)

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

- n:

  number of random values to return.

## Value

`dztbetabinom` gives the probability mass function and `rztbetabinom`
generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

By definition, this distribution only has support on the positive
integers (1, ..., n). Any zero-truncated distribution is defined as
\$\$P(X=x \| X\>0) = P(X=x) / (1 - P(X=0)),\$\$ where \\P(X=x)\\ is the
probability mass function of the corresponding untruncated distribution.

Like
[`betabinom`](https://janolefi.github.io/RTMBdist/reference/betabinom.md)
itself, this distribution provides no distribution function: the
beta-binomial cdf has no closed form and would have to be summed over
the support, which cannot be taped for automatic differentiation.

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
```
