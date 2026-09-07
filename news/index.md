# Changelog

## RTMBdist 1.1.0

- Added the hurdle (zero-altered) Poisson distribution
  ([`dhpois()`](https://janolefi.github.io/RTMBdist/reference/hpois.md),
  [`phpois()`](https://janolefi.github.io/RTMBdist/reference/hpois.md),
  [`rhpois()`](https://janolefi.github.io/RTMBdist/reference/hpois.md)),
  in which the probability of a zero is a free parameter and the
  positive counts follow the zero-truncated Poisson. Unlike
  zero-inflation, which can only add zeros, a hurdle model allows
  `zeroprob` to be smaller than the Poisson would give on its own. The
  density and distribution function are both differentiable, so
  simulation and one-step-ahead residuals via `method = "cdf"` are
  supported.

- Added the Johnson SU distribution in both the original
  parameterisation
  ([`djsu()`](https://janolefi.github.io/RTMBdist/reference/jsu.md),
  [`pjsu()`](https://janolefi.github.io/RTMBdist/reference/jsu.md),
  [`qjsu()`](https://janolefi.github.io/RTMBdist/reference/jsu.md),
  [`rjsu()`](https://janolefi.github.io/RTMBdist/reference/jsu.md)) and
  the moment parameterisation
  ([`djsu2()`](https://janolefi.github.io/RTMBdist/reference/jsu.md),
  [`pjsu2()`](https://janolefi.github.io/RTMBdist/reference/jsu.md),
  [`qjsu2()`](https://janolefi.github.io/RTMBdist/reference/jsu.md),
  [`rjsu2()`](https://janolefi.github.io/RTMBdist/reference/jsu.md)), a
  four-parameter distribution on the real line covering a wide range of
  skewness and kurtosis. In
  [`djsu2()`](https://janolefi.github.io/RTMBdist/reference/jsu.md) the
  location and scale arguments are exactly the mean and standard
  deviation. The density and distribution function are both
  differentiable, so simulation and one-step-ahead residuals are
  supported. Unlike `gamlss.dist`, the reparameterisation stays finite
  for very large `tau`, where the distribution approaches the normal.

- Added references to the primary source for each distribution derived
  from `gamlss.dist`, and cross-links between related families.
  [`pgenpois()`](https://janolefi.github.io/RTMBdist/reference/genpois.md)
  and friends previously had no references at all.

- Removed the dependency on `gamlss.dist`, which is scheduled for
  archival on CRAN. The quantile and random generation functions of the
  Box-Cox Cole-Green
  ([`qbccg()`](https://janolefi.github.io/RTMBdist/reference/bccg.md),
  [`rbccg()`](https://janolefi.github.io/RTMBdist/reference/bccg.md)),
  Box-Cox *t*
  ([`qbct()`](https://janolefi.github.io/RTMBdist/reference/bct.md),
  [`rbct()`](https://janolefi.github.io/RTMBdist/reference/bct.md)),
  Box-Cox power exponential
  ([`qbcpe()`](https://janolefi.github.io/RTMBdist/reference/bcpe.md),
  [`rbcpe()`](https://janolefi.github.io/RTMBdist/reference/bcpe.md)),
  power exponential
  ([`qpowerexp()`](https://janolefi.github.io/RTMBdist/reference/powerexp.md),
  [`rpowerexp()`](https://janolefi.github.io/RTMBdist/reference/powerexp.md),
  [`qpowerexp2()`](https://janolefi.github.io/RTMBdist/reference/powerexp.md),
  [`rpowerexp2()`](https://janolefi.github.io/RTMBdist/reference/powerexp.md)),
  Pareto
  ([`qpareto()`](https://janolefi.github.io/RTMBdist/reference/pareto.md),
  [`rpareto()`](https://janolefi.github.io/RTMBdist/reference/pareto.md)),
  generalised Poisson
  ([`pgenpois()`](https://janolefi.github.io/RTMBdist/reference/genpois.md),
  [`qgenpois()`](https://janolefi.github.io/RTMBdist/reference/genpois.md),
  [`rgenpois()`](https://janolefi.github.io/RTMBdist/reference/genpois.md))
  and exponentially modified Gaussian
  ([`qexgauss()`](https://janolefi.github.io/RTMBdist/reference/exgauss.md))
  distributions are now implemented natively. Results are unchanged
  except for the fixes below.

- Fixed argument recycling in
  [`qbccg()`](https://janolefi.github.io/RTMBdist/reference/bccg.md),
  [`qbct()`](https://janolefi.github.io/RTMBdist/reference/bct.md),
  [`qbcpe()`](https://janolefi.github.io/RTMBdist/reference/bcpe.md),
  [`qgenpois()`](https://janolefi.github.io/RTMBdist/reference/genpois.md)
  and
  [`pgenpois()`](https://janolefi.github.io/RTMBdist/reference/genpois.md).
  Evaluating a single quantile against vectorised parameters previously
  collapsed the result to length one, and a parameter vector shorter
  than `x` silently truncated it. These functions now return one value
  per recycled argument tuple.

- Fixed [`qbct()`](https://janolefi.github.io/RTMBdist/reference/bct.md)
  and [`rbct()`](https://janolefi.github.io/RTMBdist/reference/bct.md),
  which failed with an error when `mu`, `sigma` or `tau` was supplied as
  a vector.

- [`qpareto()`](https://janolefi.github.io/RTMBdist/reference/pareto.md)
  now honours `lower.tail`, which was previously accepted but ignored.

- `log.p = TRUE` now works in
  [`qpareto()`](https://janolefi.github.io/RTMBdist/reference/pareto.md),
  [`qbct()`](https://janolefi.github.io/RTMBdist/reference/bct.md) and
  [`qgenpois()`](https://janolefi.github.io/RTMBdist/reference/genpois.md).
  These previously validated `p` before transforming it back from the
  log scale, so the argument could not be used.

- [`qexgauss()`](https://janolefi.github.io/RTMBdist/reference/exgauss.md)
  is substantially more accurate. It inverts
  [`pexgauss()`](https://janolefi.github.io/RTMBdist/reference/exgauss.md)
  with a tighter convergence tolerance, reducing the round-trip error
  `|p(q(p)) - p|` from roughly 1e-6 to roughly 1e-14.

- [`pgenpois()`](https://janolefi.github.io/RTMBdist/reference/genpois.md)
  no longer falls back to the Poisson distribution for `phi < 1e-4` and
  evaluates the generalised Poisson distribution function across the
  whole parameter range.

- [`pbetaprime()`](https://janolefi.github.io/RTMBdist/reference/betaprime.md),
  [`pinvchisq()`](https://janolefi.github.io/RTMBdist/reference/invchisq.md)
  and
  [`pinvgamma()`](https://janolefi.github.io/RTMBdist/reference/invgamma.md)
  are now differentiable, so one-step-ahead residuals via
  `method = "cdf"` are available for the beta prime, inverse chi-squared
  and inverse gamma distributions.

- Added the Bell distribution
  ([`dbell()`](https://janolefi.github.io/RTMBdist/reference/bell.md),
  [`pbell()`](https://janolefi.github.io/RTMBdist/reference/bell.md),
  [`qbell()`](https://janolefi.github.io/RTMBdist/reference/bell.md),
  [`rbell()`](https://janolefi.github.io/RTMBdist/reference/bell.md))
  and its mean parameterisation
  ([`dbell2()`](https://janolefi.github.io/RTMBdist/reference/bell2.md),
  [`pbell2()`](https://janolefi.github.io/RTMBdist/reference/bell2.md),
  [`qbell2()`](https://janolefi.github.io/RTMBdist/reference/bell2.md),
  [`rbell2()`](https://janolefi.github.io/RTMBdist/reference/bell2.md)),
  a one-parameter distribution for overdispersed counts. Log Bell
  numbers are used throughout and cached, so the density stays finite
  well past the point at which the Bell numbers themselves overflow
  double precision. One-step-ahead residuals are supported via
  `method = "cdf"`.

- Added
  [`lambertW()`](https://janolefi.github.io/RTMBdist/reference/lambertW.md),
  an AD-compatible implementation of the principal branch of the Lambert
  W function. Derivatives of all orders are exact, so it can be used in
  models fitted by Laplace approximation.

- Added the Conway-Maxwell-binomial distribution
  ([`dcombinom()`](https://janolefi.github.io/RTMBdist/reference/combinom.md),
  [`pcombinom()`](https://janolefi.github.io/RTMBdist/reference/combinom.md),
  [`qcombinom()`](https://janolefi.github.io/RTMBdist/reference/combinom.md),
  [`rcombinom()`](https://janolefi.github.io/RTMBdist/reference/combinom.md)),
  a binomial generalisation with an additional dispersion parameter
  `nu`. The density and distribution function are differentiable with
  respect to `x`, so one-step-ahead residuals are supported.

## RTMBdist 0.1.0

CRAN release: 2025-10-07

- First release on CRAN.
