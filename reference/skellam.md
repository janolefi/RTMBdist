# Skellam distribution

Probability mass function, distribution function, and random generation
for the Skellam distribution.

## Usage

``` r
dskellam(x, mu1, mu2, log = FALSE)

pskellam(q, mu1, mu2, lower.tail = TRUE, log.p = FALSE)

rskellam(n, mu1, mu2)
```

## Arguments

- x:

  integer vector of counts

- mu1, mu2:

  Poisson means

- log:

  logical; return log-density if TRUE

- q:

  vector of quantiles

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le q\]\\,
  otherwise \\P\[X \> q\]\\.

- log.p:

  logical; if `TRUE`, probabilities are returned on the log scale.

- n:

  number of random values to return.

## Value

`dskellam` gives the probability mass function, `pskellam` gives the
distribution function, and `rskellam` generates random deviates.

## Details

The Skellam distribution is the distribution of the difference of two
Poisson random variables. Specifically, if \\X_1 \sim
\text{Pois}(\mu_1)\\ and \\X_2 \sim \text{Pois}(\mu_2)\\, then \\X_1 -
X_2 \sim \text{Skellam}(\mu_1, \mu_2)\\.

This implementation of `dskellam` allows for automatic differentiation
with `RTMB`.

The distribution function has no closed form, and as the support is
unbounded in both directions, summing the probability mass function
would need a range that depends on the parameters. Instead, `pskellam`
uses the identity \$\$P(X \le q;\\ \mu_1, \mu_2) = \int\_{\mu_1}^\infty
P(X = q;\\ t, \mu_2)\\ dt,\$\$ which follows from the relation between
the Poisson and gamma distribution functions, and integrates it
numerically with the AD-compatible
[`integrate`](https://rdrr.io/pkg/RTMB/man/ADintegrate.html) of `RTMB`.
By symmetry, the upper tail is \\P(X \> q;\\ \mu_1, \mu_2) = P(X \le
-q - 1;\\ \mu_2, \mu_1)\\, and whichever tail is smaller is integrated
directly, which keeps small tail probabilities accurate. `pskellam` is
AD-compatible in `mu1` and `mu2`, while `q` must be numeric data. This
is also what one-step-ahead (OSA) residuals via
`RTMB::`[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
need, so these are supported, e.g. with `method = "cdf"` and
`discrete = TRUE`.

## Examples

``` r
x <- rskellam(1, 2, 3)
d <- dskellam(x, 2, 3)
p <- pskellam(x, 2, 3)
```
