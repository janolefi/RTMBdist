# wrapped Cauchy distribution

Density, distribution function, quantile function, and random generation
for the wrapped Cauchy distribution.

## Usage

``` r
dwrpcauchy(x, mu = 0, rho, log = FALSE)

pwrpcauchy(q, mu = 0, rho, from = NULL, lower.tail = TRUE, log.p = FALSE)

qwrpcauchy(p, mu = 0, rho, from = NULL, lower.tail = TRUE, log.p = FALSE)

rwrpcauchy(n, mu = 0, rho, wrap = TRUE)
```

## Arguments

- x, q:

  vector of angles measured in radians at which to evaluate the density
  or distribution function.

- mu:

  mean direction of the distribution measured in radians.

- rho:

  concentration parameter of the distribution, must be in the interval
  from 0 to 1.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- from:

  origin, in radians, at which the circle is cut open for the
  distribution and quantile functions. If `NULL` (default), it is set to
  `mu - pi`.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities.

- n:

  number of random values to return.

- wrap:

  logical; if `TRUE`, generated angles are wrapped to the interval from
  -pi to pi, otherwise they lie in the interval from `mu - pi` to
  `mu + pi`.

## Value

`dwrpcauchy` gives the density, `pwrpcauchy` gives the distribution
function, `qwrpcauchy` gives the quantile function, and `rwrpcauchy`
generates random deviates.

## Details

`dwrpcauchy` and `pwrpcauchy` allow for automatic differentiation with
`RTMB`.

\$\$f(x;\\\mu,\rho) = \frac{1}{2\pi}\cdot\frac{1-\rho^2}{1 + \rho^2 -
2\rho\cos(x-\mu)}.\$\$

A circular distribution has no smallest angle, so its distribution
function depends on where the circle is cut open. `pwrpcauchy` cuts it
at the antipode of the mean direction, \\\mu - \pi\\, and is then
available in closed form: \$\$F(q) = P(\mu - \pi \< X \le q) =
\frac{1}{2} +
\frac{1}{\pi}\arctan\left(\frac{1+\rho}{1-\rho}\tan\frac{q-\mu}{2}\right),
\quad \mu - \pi \< q \le \mu + \pi,\$\$ so that \\F(\mu) = 1/2\\. This
is the default, `from = NULL`, and the same origin as the default of
[`pvm`](https://janolefi.github.io/RTMBdist/reference/vm.md).

Angles outside \\(\mu - \pi, \mu + \pi\]\\ are wrapped onto this
interval first, so `pwrpcauchy` is \\2\pi\\-periodic in `q`. As a
consequence, the default origin moves with `mu`: for angles on \\\[-\pi,
\pi\]\\ and \\\mu \neq 0\\, `pwrpcauchy` is not monotone over that range
but drops from 1 back to 0 at \\\mu \pm \pi\\.

A different origin is set with `from`, giving \\P(\mathrm{from} \< X \le
q) = (F(q) - F(\mathrm{from})) \bmod 1\\ for \\\mathrm{from} \< q \le
\mathrm{from} + 2\pi\\. A fixed `from` is needed whenever the
distribution function is averaged over different values of `mu`, for
example over the states of a hidden Markov model or over a random
effect: with the default, each value of `mu` cuts the circle at a
different place, and the average of these distribution functions is not
the distribution function of the mixture.

`qwrpcauchy` is the inverse of `pwrpcauchy` and returns angles in
\\\[\mathrm{from}, \mathrm{from} + 2\pi\]\\, by default \\\[\mu - \pi,
\mu + \pi\]\\, not wrapped to \\\[-\pi, \pi\]\\. The latter also holds
for `rwrpcauchy` with `wrap = FALSE`.

## See also

[vm](https://janolefi.github.io/RTMBdist/reference/vm.md);
[`cjw()`](https://janolefi.github.io/RTMBdist/reference/cjw.md) and
[`cfold()`](https://janolefi.github.io/RTMBdist/reference/cfold.md) for
circular-linear copulas joining turning angles and step lengths.

## Examples

``` r
set.seed(1)
x <- rwrpcauchy(10, 0, 0.5)
d <- dwrpcauchy(x, 0, 0.5)
p <- pwrpcauchy(x, 0, 0.5)
q <- qwrpcauchy(p, 0, 0.5)
```
