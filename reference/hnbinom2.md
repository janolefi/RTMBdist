# Reparameterised hurdle negative binomial distribution

Probability mass function, distribution function, and random generation
for the hurdle (zero-altered) negative binomial distribution,
parameterised by the mean of the untruncated negative binomial.

## Usage

``` r
dhnbinom2(x, mu, size, zeroprob = 0.5, log = FALSE)

phnbinom2(q, mu, size, zeroprob = 0.5, lower.tail = TRUE, log.p = FALSE)

rhnbinom2(n, mu, size, zeroprob = 0.5)
```

## Arguments

- x, q:

  integer vector of counts

- mu:

  mean of the untruncated negative binomial, must be strictly positive

- size:

  dispersion parameter, must be strictly positive

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

`dhnbinom2` gives the probability mass function, `phnbinom2` gives the
distribution function, and `rhnbinom2` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

This is
[`hnbinom`](https://janolefi.github.io/RTMBdist/reference/hnbinom.md)
with the success probability replaced by \$\$\pi =
\frac{\mathrm{size}}{\mathrm{size} + \mu},\$\$ so that \\\mu\\ is the
mean of the *untruncated* negative binomial, whose variance is \\\mu +
\mu^2/\mathrm{size}\\. Note that \\\mu\\ is not the mean of the hurdle
distribution itself, which also depends on `zeroprob`.

As for all hurdle distributions, `zeroprob` is exactly the probability
of observing a zero and may be larger *or* smaller than the negative
binomial would give on its own.

## References

Mullahy, J. (1986) Specification and testing of some modified count data
models. Journal of Econometrics, 33, 341-365.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[hnbinom](https://janolefi.github.io/RTMBdist/reference/hnbinom.md),
[hpois](https://janolefi.github.io/RTMBdist/reference/hpois.md),
[nbinom2](https://janolefi.github.io/RTMBdist/reference/nbinom2.md),
[zinbinom2](https://janolefi.github.io/RTMBdist/reference/zinbinom2.md),
[ztnbinom2](https://janolefi.github.io/RTMBdist/reference/ztnbinom2.md)

## Examples

``` r
set.seed(123)
x <- rhnbinom2(5, mu = 3, size = 2, zeroprob = 0.3)
d <- dhnbinom2(x, mu = 3, size = 2, zeroprob = 0.3)
p <- phnbinom2(x, mu = 3, size = 2, zeroprob = 0.3)
```
