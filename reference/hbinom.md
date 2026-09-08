# Hurdle binomial distribution

Probability mass function, distribution function, and random generation
for the hurdle (zero-altered) binomial distribution.

## Usage

``` r
dhbinom(x, size, prob, zeroprob = 0.5, log = FALSE)

phbinom(q, size, prob, zeroprob = 0.5, lower.tail = TRUE, log.p = FALSE)

rhbinom(n, size, prob, zeroprob = 0.5)
```

## Arguments

- x, q:

  integer vector of counts

- size:

  number of trials (zero or more)

- prob:

  probability of success on each trial

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

`dhbinom` gives the probability mass function, `phbinom` gives the
distribution function, and `rhbinom` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

A hurdle distribution models the zeros and the positive counts as two
separate processes: the probability of a zero is a free parameter, and
the positive counts follow the corresponding zero-truncated
distribution. Writing \\p_0\\ for `zeroprob`, \$\$P(X = 0) = p_0, \qquad
P(X = x) = (1 - p_0)\\\frac{P\_{\mathrm{Bin}}(x;\\n,\pi)}{1 - \pi_0},
\quad x = 1, \ldots, n.\$\$ where \\\pi_0 =
P\_{\mathrm{Bin}}(0;\\n,\pi)\\ is the probability of a zero under the
ordinary binomial.

Unlike zero-inflation, which can only add zeros to those the binomial
already produces, `zeroprob` here is exactly the probability of a zero
and may be larger *or* smaller than \\\pi_0\\. The two coincide with the
ordinary binomial when `zeroprob` equals \\\pi_0\\.

## References

Mullahy, J. (1986) Specification and testing of some modified count data
models. Journal of Econometrics, 33, 341-365.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[hpois](https://janolefi.github.io/RTMBdist/reference/hpois.md),
[hnbinom](https://janolefi.github.io/RTMBdist/reference/hnbinom.md),
[zibinom](https://janolefi.github.io/RTMBdist/reference/zibinom.md),
[ztbinom](https://janolefi.github.io/RTMBdist/reference/ztbinom.md)

## Examples

``` r
set.seed(123)
x <- rhbinom(5, size = 10, prob = 0.3, zeroprob = 0.4)
d <- dhbinom(x, size = 10, prob = 0.3, zeroprob = 0.4)
p <- phbinom(x, size = 10, prob = 0.3, zeroprob = 0.4)
```
