# Hurdle beta-binomial distribution

Probability mass function and random generation for the hurdle
(zero-altered) beta-binomial distribution.

## Usage

``` r
dhbetabinom(x, size, shape1, shape2, zeroprob = 0.5, log = FALSE)

rhbetabinom(n, size, shape1, shape2, zeroprob = 0.5)
```

## Arguments

- x:

  integer vector of counts

- size:

  number of trials (zero or more)

- shape1, shape2:

  positive shape parameters of the mixing beta distribution

- zeroprob:

  probability of a zero, between 0 and 1

- log:

  logical; return log-density if TRUE

- n:

  number of random values to return.

## Value

`dhbetabinom` gives the probability mass function and `rhbetabinom`
generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

A hurdle distribution models the zeros and the positive counts as two
separate processes: the probability of a zero is a free parameter, and
the positive counts follow the corresponding zero-truncated
distribution. Writing \\p_0\\ for `zeroprob`, \$\$P(X = 0) = p_0, \qquad
P(X = x) = (1 - p_0)\\\frac{P\_{\mathrm{BB}}(x;\\n,a,b)}{1 - \pi_0},
\quad x = 1, \ldots, n.\$\$ where \\\pi_0 =
P\_{\mathrm{BB}}(0;\\n,a,b)\\ is the probability of a zero under the
ordinary beta-binomial.

Unlike zero-inflation, which can only add zeros, `zeroprob` here is
exactly the probability of observing a zero and may be larger *or*
smaller than the beta-binomial would give on its own.

Like
[`betabinom`](https://janolefi.github.io/RTMBdist/reference/betabinom.md)
itself, this distribution provides no distribution function: the
beta-binomial cdf has no closed form and would have to be summed over
the support, which cannot be taped for automatic differentiation.

## See also

[betabinom](https://janolefi.github.io/RTMBdist/reference/betabinom.md),
[zibetabinom](https://janolefi.github.io/RTMBdist/reference/zibetabinom.md),
[ztbetabinom](https://janolefi.github.io/RTMBdist/reference/ztbetabinom.md),
[hbinom](https://janolefi.github.io/RTMBdist/reference/hbinom.md)

## Examples

``` r
set.seed(123)
x <- rhbetabinom(5, size = 10, shape1 = 2, shape2 = 3, zeroprob = 0.4)
d <- dhbetabinom(x, size = 10, shape1 = 2, shape2 = 3, zeroprob = 0.4)
```
