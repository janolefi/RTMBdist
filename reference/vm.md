# von Mises distribution

Density, distribution function, and random generation for the von Mises
distribution.

## Usage

``` r
dvm(x, mu = 0, kappa = 1, log = FALSE)

pvm(
  q,
  mu = 0,
  kappa = 1,
  from = NULL,
  tol = 1e-20,
  lower.tail = TRUE,
  log.p = FALSE
)

rvm(n, mu = 0, kappa = 1, wrap = TRUE)
```

## Arguments

- x, q:

  vector of angles measured in radians at which to evaluate the density
  function.

- mu:

  mean direction of the distribution measured in radians.

- kappa:

  non-negative numeric value for the concentration parameter of the
  distribution.

- log:

  logical; if `TRUE`, densities are returned on the log scale.

- from:

  value from which the integration for CDF starts. If `NULL`, is set to
  `mu - pi`.

- tol:

  the precision in evaluating the distribution function, ignored in AD
  context.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- log.p:

  logical; if `TRUE`, probabilities are returned as \\\log(p)\\.

- n:

  number of random values to return.

- wrap:

  logical; if `TRUE`, generated angles are wrapped to the interval from
  -pi to pi.

## Value

`dvm` gives the density, `pvm` gives the distribution function, and
`rvm` generates random deviates.

## Details

This implementation of `dvm` allows for automatic differentiation with
`RTMB`. `rvm` and `pvm` are simply wrappers of the corresponding
functions from `circular`. When called during AD taping, `pvm` instead
integrates the density numerically with the AD-compatible
[`integrate`](https://rdrr.io/pkg/RTMB/man/ADintegrate.html) of `RTMB`,
which makes it AD-compatible in `q`, `mu`, `kappa` and `from`. The `tol`
argument is then ignored.

\$\$f(x;\\\mu,\kappa) = \frac{\exp(\kappa\cos(x-\mu))}{2\pi\\
I_0(\kappa)},\$\$ where \\I_0\\ is the modified Bessel function of the
first kind of order 0.

A circular distribution has no smallest angle, so its distribution
function depends on where the circle is cut open. By default, `pvm` cuts
it at the antipode of the mean direction, \\\mu - \pi\\, so that
\\F(\mu) = 1/2\\. A different origin is set with `from`. A fixed `from`
is needed whenever the distribution function is averaged over different
values of `mu`, for example over the states of a hidden Markov model or
over a random effect: with the default, each value of `mu` cuts the
circle at a different place, and the average of these distribution
functions is not the distribution function of the mixture.

**OSA residuals:** `dvm` supports one-step-ahead (OSA) quantile
residuals via
`RTMB::`[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html).
For the methods based on the distribution function, such as
`method = "cdf"`, the circle is cut at the fixed origin \\-\pi\\, i.e.
the residuals are based on `pvm(x, mu, kappa, from = -pi)` rather than
the default origin \\\mu - \pi\\. OSA residuals are computed from the
predictive distribution function, which averages the distribution
function over hidden states or random effects, and this is only valid
with an origin that does not depend on `mu` (see above). Hence the
residuals are valid for all models, but their interpretation depends on
the data: for turning angles, with `mu` close to 0, the cut at
\\\pm\pi\\ corresponds to a reversal and the residuals increase with the
turning angle. For directions with `mu` far from 0, angles close to
\\\pm\pi\\ can give large residuals of either sign, even when they are
close to the mean direction. For `method = "oneStepGeneric"`, set
`range = c(-pi, pi)` in
[`oneStepPredict()`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html),
so that the density is integrated from the same origin.

## See also

[wrpcauchy](https://janolefi.github.io/RTMBdist/reference/wrpcauchy.md);
[`cjw()`](https://janolefi.github.io/RTMBdist/reference/cjw.md) and
[`cfold()`](https://janolefi.github.io/RTMBdist/reference/cfold.md) for
circular-linear copulas joining turning angles and step lengths.

## Examples

``` r
set.seed(1)
x <- rvm(10, 0, 1)
d <- dvm(x, 0, 1)
p <- pvm(x, 0, 1)
```
