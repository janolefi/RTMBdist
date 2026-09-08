# Hurdle negative binomial distribution

Probability mass function, distribution function, and random generation
for the hurdle (zero-altered) negative binomial distribution.

## Usage

``` r
dhnbinom(x, size, prob, zeroprob = 0.5, log = FALSE)

phnbinom(q, size, prob, zeroprob = 0.5, lower.tail = TRUE, log.p = FALSE)

rhnbinom(n, size, prob, zeroprob = 0.5)
```

## Arguments

- x, q:

  integer vector of counts

- size:

  dispersion parameter, must be strictly positive

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

`dhnbinom` gives the probability mass function, `phnbinom` gives the
distribution function, and `rhnbinom` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

A hurdle distribution models the zeros and the positive counts as two
separate processes: the probability of a zero is a free parameter, and
the positive counts follow the corresponding zero-truncated
distribution. Writing \\p_0\\ for `zeroprob`, \$\$P(X = 0) = p_0, \qquad
P(X = x) = (1 - p_0)\\\frac{P\_{\mathrm{NB}}(x;\\r,\pi)}{1 - \pi_0},
\quad x = 1, 2, \ldots\$\$ where \\\pi_0 = P\_{\mathrm{NB}}(0;\\r,\pi)\\
is the probability of a zero under the ordinary negative binomial.

Unlike zero-inflation, which can only add zeros to those the negative
binomial already produces, `zeroprob` here is exactly the probability of
a zero and may be larger *or* smaller than \\\pi_0\\. The two coincide
with the ordinary negative binomial when `zeroprob` equals \\\pi_0\\.

## References

Mullahy, J. (1986) Specification and testing of some modified count data
models. Journal of Econometrics, 33, 341-365.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[hnbinom2](https://janolefi.github.io/RTMBdist/reference/hnbinom2.md),
[hpois](https://janolefi.github.io/RTMBdist/reference/hpois.md),
[hbinom](https://janolefi.github.io/RTMBdist/reference/hbinom.md),
[zinbinom](https://janolefi.github.io/RTMBdist/reference/zinbinom.md),
[ztnbinom](https://janolefi.github.io/RTMBdist/reference/ztnbinom.md)

## Examples

``` r
set.seed(123)
x <- rhnbinom(5, size = 2, prob = 0.4, zeroprob = 0.3)
d <- dhnbinom(x, size = 2, prob = 0.4, zeroprob = 0.3)
p <- phnbinom(x, size = 2, prob = 0.4, zeroprob = 0.3)
```
