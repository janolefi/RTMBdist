# Bell distribution

Probability mass function, distribution function, quantile function, and
random generation for the Bell distribution.

## Usage

``` r
dbell(x, theta, log = FALSE)

pbell(q, theta, lower.tail = TRUE, log.p = FALSE)

qbell(p, theta, lower.tail = TRUE, log.p = FALSE)

rbell(n, theta)
```

## Arguments

- x, q:

  integer vector of counts

- theta:

  vector of positive Bell parameters

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return.

## Value

`dbell` gives the probability mass function, `pbell` gives the
distribution function, `qbell` gives the quantile function, and `rbell`
generates random deviates.

## Details

This implementation of `dbell` and `pbell` allows for automatic
differentiation with `RTMB` with respect to `theta`.

The Bell distribution (Castellares et al. 2018) is a one-parameter
distribution for overdispersed counts with probability mass function

\$\$P(X = x;\\ \theta) = \frac{e^{1 - e^{\theta}}\\ \theta^{x}\\
B_x}{x!}, \quad x = 0, 1, 2, \ldots,\$\$

for \\\theta \> 0\\, where \\B_x\\ is the \\x\\-th Bell number, i.e. the
number of ways a set of \\x\\ elements can be partitioned into non-empty
subsets.

Its mean and variance are \$\$E(X) = \theta e^{\theta}, \qquad
\mathrm{Var}(X) = \theta e^{\theta} (1 + \theta),\$\$ so the dispersion
index is \\\mathrm{Var}(X) / E(X) = 1 + \theta \> 1\\ and the
distribution is always overdispersed. Note that, unlike the negative
binomial or generalised Poisson distributions, the Bell distribution has
no separate dispersion parameter: the amount of overdispersion is tied
to the mean. As \\\theta \to 0\\ it approaches the Poisson distribution.

The distribution arises as a compound Poisson sum \$\$X =
\sum\_{i=1}^{N} Y_i, \qquad N \sim \mathrm{Pois}(e^{\theta} - 1), \quad
Y_i \sim \mathrm{ztPois}(\theta),\$\$ with the \\Y_i\\ independent of
\\N\\, and `rbell` uses exactly this representation. It is therefore
infinitely divisible, and a member of the one-parameter exponential
family with natural parameter \\\log\theta\\ and sufficient statistic
\\x\\.

The Bell numbers grow faster than any exponential and overflow double
precision at \\x = 219\\, so \\\log B_x\\ is used throughout rather than
\\B_x\\. As `x` is data, \\\log B_x\\ is constant with respect to the
parameters and is cached across calls.

Neither the distribution function nor the quantile function has a closed
form; both are obtained by summing the probability mass function, with
`qbell` choosing its summation range automatically from the mean and
variance.

See [`bell2`](https://janolefi.github.io/RTMBdist/reference/bell2.md)
for the parameterisation by the mean.

## References

Castellares, F., Ferrari, S. L. P., and Lemonte, A. J. (2018). On the
Bell distribution and its associated regression model for count data.
*Applied Mathematical Modelling* 56, 172-185.
doi:10.1016/j.apm.2017.12.014

## Examples

``` r
set.seed(123)
x <- rbell(1, 1)
d <- dbell(x, 1)
p <- pbell(x, 1)
q <- qbell(p, 1)

# mean and variance
theta <- 0.8
xs <- 0:100
sum(xs * dbell(xs, theta)) # theta * exp(theta)
#> [1] 1.780433
```
