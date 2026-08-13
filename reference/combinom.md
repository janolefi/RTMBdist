# Conway-Maxwell-binomial distribution

Probability mass function, distribution function, quantile function, and
random generation for the Conway-Maxwell-binomial (CMB) distribution.

## Usage

``` r
dcombinom(x, size, prob, nu = 1, log = FALSE)

pcombinom(q, size, prob, nu = 1, lower.tail = TRUE, log.p = FALSE)

qcombinom(p, size, prob, nu = 1, lower.tail = TRUE, log.p = FALSE)

rcombinom(n, size, prob, nu = 1)
```

## Arguments

- x, q:

  integer vector of counts in \\\\0, 1, \ldots, \\`size`\\\\\\

- size:

  vector of numbers of trials (non-negative integers)

- prob:

  vector of success probabilities in \\(0, 1)\\

- nu:

  vector of dispersion parameters; `nu = 1` gives the binomial
  distribution, `nu > 1` under-dispersion and `nu < 1` over-dispersion.
  May be negative.

- log, log.p:

  logical; if `TRUE`, probabilities/densities are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dcombinom` gives the probability mass function, `pcombinom` gives the
distribution function, `qcombinom` gives the quantile function, and
`rcombinom` generates random deviates.

## Details

This implementation of `dcombinom` and `pcombinom` allows for automatic
differentiation with `RTMB`, including differentiation with respect to
`x`, such that one-step-ahead (OSA) residuals are supported.

The CMB distribution generalises the binomial distribution by an
additional dispersion parameter \\\nu\\ in the same way that the
Conway-Maxwell-Poisson distribution generalises the Poisson
distribution. Its probability mass function is

\$\$P(X = x;\\ n, p, \nu) = \frac{1}{Z(n, p, \nu)} \binom{n}{x}^{\nu}
p^x (1-p)^{n-x}, \quad x = 0, 1, \ldots, n,\$\$

with normalising constant

\$\$Z(n, p, \nu) = \sum\_{k=0}^{n} \binom{n}{k}^{\nu} p^k
(1-p)^{n-k}.\$\$

For \\\nu = 1\\ this reduces to the binomial distribution. Values \\\nu
\> 1\\ give under-dispersion and \\\nu \< 1\\ over-dispersion relative
to a binomial distribution with the same mean. As the support is finite,
\\Z\\ converges for every real \\\nu\\, so \\\nu\\ is not restricted to
be positive; negative values yield strongly over-dispersed, U-shaped
distributions.

The distribution arises as the sum of \\n\\ exchangeable, possibly
associated Bernoulli variables, where \\\nu\\ controls the association.
Note that `prob` is *not* the mean divided by `size` unless \\\nu = 1\\;
the mean has no closed form for general \\\nu\\ and must be obtained by
summation over the support.

## References

Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S., and Boatwright, P.
(2005). A useful distribution for fitting discrete data: revival of the
Conway-Maxwell-Poisson distribution. *Journal of the Royal Statistical
Society: Series C* 54(1), 127-142.

Kadane, J. B. (2016). Sums of possibly associated Bernoulli variables:
the Conway-Maxwell-binomial distribution. *Bayesian Analysis* 11(2),
403-420.

<https://en.wikipedia.org/wiki/Conway-Maxwell-binomial_distribution>

## Examples

``` r
set.seed(123)
x <- rcombinom(1, size = 10, prob = 0.4, nu = 1.5)
d <- dcombinom(x, 10, 0.4, 1.5)
p <- pcombinom(x, 10, 0.4, 1.5)
q <- qcombinom(p, 10, 0.4, 1.5)

# nu = 1 recovers the binomial distribution
all.equal(dcombinom(0:10, 10, 0.3, 1), dbinom(0:10, 10, 0.3))
#> [1] TRUE
```
