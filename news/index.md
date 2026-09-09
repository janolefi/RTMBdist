# Changelog

## RTMBdist 1.2.0

- Added the half-t
  ([`dhalft()`](https://janolefi.github.io/RTMBdist/reference/halft.md))
  and half-Cauchy
  ([`dhalfcauchy()`](https://janolefi.github.io/RTMBdist/reference/halfcauchy.md))
  distributions, each with matching `p`, `q` and `r` functions. These
  are the standard weakly informative priors for the standard deviation
  of a hierarchical model, so they are aimed squarely at models fitted
  by the Laplace approximation. Both distribution functions are
  differentiable, and
  [`phalft()`](https://janolefi.github.io/RTMBdist/reference/halft.md)
  differentiates with respect to `df` as well as `sigma`, so the degrees
  of freedom can be estimated rather than fixed; one-step-ahead
  residuals via `method = "cdf"` are supported. The half-Cauchy is the
  half-t with `df = 1`, and the half-normal is
  [`dfoldnorm()`](https://janolefi.github.io/RTMBdist/reference/foldnorm.md)
  with `mu = 0`, which the half-t approaches as `df` grows.

- Added the Yule-Simon
  ([`dyules()`](https://janolefi.github.io/RTMBdist/reference/yules.md))
  and Waring
  ([`dwaring()`](https://janolefi.github.io/RTMBdist/reference/waring.md))
  distributions, two classical long-tailed count laws, each with
  matching `p` and `r` functions. Both are beta-geometric special cases
  of the beta-negative binomial, with `size` fixed at one, and that
  restriction is what gives them a closed-form distribution function
  where the general beta-negative binomial has none; one-step-ahead
  residuals via `method = "cdf"` are therefore available for these two.
  [`dyules()`](https://janolefi.github.io/RTMBdist/reference/yules.md)
  follows `VGAM` and is supported on the positive integers, while
  [`dwaring()`](https://janolefi.github.io/RTMBdist/reference/waring.md)
  follows the `WARING` family of `gamlss.dist`, starts at zero and has
  `mu` as its mean. The `YULE` family of `gamlss.dist` is
  `dwaring(x, mu, mu)`.

- Added the beta-negative binomial distribution
  ([`dbnbinom()`](https://janolefi.github.io/RTMBdist/reference/bnbinom.md))
  and its mean parameterisation
  ([`dbnbinom2()`](https://janolefi.github.io/RTMBdist/reference/bnbinom2.md)),
  each with a matching `r` function. It is to the negative binomial what
  the beta-binomial is to the binomial, and its extra beta layer gives a
  considerably heavier tail. In
  [`dbnbinom2()`](https://janolefi.github.io/RTMBdist/reference/bnbinom2.md),
  which follows the `BNB` family of `gamlss.dist`, `mu` is exactly the
  mean; this is the more stable parameterisation to estimate in, because
  the original one has a long likelihood ridge along which `size` and
  `shape2` trade off. Like the beta-binomial, neither has a distribution
  function, so one-step-ahead residuals are not available.

- Added the extreme value distributions: the generalised extreme value
  distribution
  ([`dgev()`](https://janolefi.github.io/RTMBdist/reference/gev.md)),
  the generalised Pareto distribution
  ([`dgpd()`](https://janolefi.github.io/RTMBdist/reference/gpd.md)) and
  the Frechet distribution
  ([`dfrechet()`](https://janolefi.github.io/RTMBdist/reference/frechet.md)),
  each with matching `p`, `q` and `r` functions. The first two cover
  their three shape regimes with a single expression rather than a
  branch on the sign of `xi`, so the derivative with respect to the
  shape is exact at `xi = 0`, which is the usual starting value when the
  shape is estimated. Densities and distribution functions are
  differentiable, so simulation and one-step-ahead residuals via
  `method = "cdf"` are supported. Unlike `VGAM`, `evd` and `extraDistr`,
  [`dgpd()`](https://janolefi.github.io/RTMBdist/reference/gpd.md)
  returns `1 / sigma` rather than zero at the threshold itself, matching
  [`stats::dexp()`](https://rdrr.io/r/stats/Exponential.html) at zero.

- `RTMBdist` no longer masks anything in `stats`. The AD-compatible
  replacements for [`stats::pt()`](https://rdrr.io/r/stats/TDist.html),
  [`stats::plnorm()`](https://rdrr.io/r/stats/Lognormal.html),
  [`stats::dgeom()`](https://rdrr.io/r/stats/Geometric.html) and
  [`stats::pgeom()`](https://rdrr.io/r/stats/Geometric.html) are
  exported as
  [`pt.ad()`](https://janolefi.github.io/RTMBdist/reference/t2.md),
  [`plnorm.ad()`](https://janolefi.github.io/RTMBdist/reference/zilnorm.md),
  [`dgeom.ad()`](https://janolefi.github.io/RTMBdist/reference/geom.ad.md)
  and
  [`pgeom.ad()`](https://janolefi.github.io/RTMBdist/reference/geom.ad.md),
  and are reached through internal S4 generics that dispatch on the
  argument classes: plain numeric input goes to the `stats` versions, AD
  variables to the `.ad` versions. Previously
  [`pt()`](https://rdrr.io/r/stats/TDist.html) was exported with a
  reduced argument list, so `pt(q, df, lower.tail = FALSE)` failed for
  anyone who had loaded the package. As a side effect `lower.tail` and
  `log.p` now also work for [`pt()`](https://rdrr.io/r/stats/TDist.html)
  under automatic differentiation.

- Added the zero-inflated
  ([`dzigeom()`](https://janolefi.github.io/RTMBdist/reference/zigeom.md)),
  zero-truncated
  ([`dztgeom()`](https://janolefi.github.io/RTMBdist/reference/ztgeom.md))
  and hurdle
  ([`dhgeom()`](https://janolefi.github.io/RTMBdist/reference/hgeom.md))
  geometric distributions. These build on `RTMB`‘s AD-compatible
  negative binomial with `size = 1`, so `stats`’ own
  [`dgeom()`](https://rdrr.io/r/stats/Geometric.html) and friends are
  left untouched.

- Added the zero-inflated, zero-truncated and hurdle beta-binomial
  distributions
  ([`dzibetabinom()`](https://janolefi.github.io/RTMBdist/reference/zibetabinom.md),
  [`dztbetabinom()`](https://janolefi.github.io/RTMBdist/reference/ztbetabinom.md),
  [`dhbetabinom()`](https://janolefi.github.io/RTMBdist/reference/hbetabinom.md)).
  Like
  [`dbetabinom()`](https://janolefi.github.io/RTMBdist/reference/betabinom.md)
  itself these have no distribution function, since the beta-binomial
  cdf has no closed form.

- The documentation of the continuous zero-inflated distributions now
  explains that, because the continuous part places no mass at zero,
  `zeroprob` is exactly the probability of a zero and zero-inflation
  coincides with a hurdle model; GAMLSS calls these zero-adjusted rather
  than zero-inflated.

- Added the hurdle (zero-altered) count distributions: Poisson
  ([`dhpois()`](https://janolefi.github.io/RTMBdist/reference/hpois.md)),
  binomial
  ([`dhbinom()`](https://janolefi.github.io/RTMBdist/reference/hbinom.md)),
  negative binomial
  ([`dhnbinom()`](https://janolefi.github.io/RTMBdist/reference/hnbinom.md))
  and its mean parameterisation
  ([`dhnbinom2()`](https://janolefi.github.io/RTMBdist/reference/hnbinom2.md)),
  each with matching `p` and `r` functions. In a hurdle distribution the
  probability of a zero is a free parameter and the positive counts
  follow the corresponding zero-truncated distribution. Unlike
  zero-inflation, which can only add zeros, a hurdle model allows
  `zeroprob` to be smaller than the Poisson would give on its own. The
  density and distribution function are both differentiable, so
  simulation and one-step-ahead residuals via `method = "cdf"` are
  supported.

- Fixed
  [`pztnbinom()`](https://janolefi.github.io/RTMBdist/reference/ztnbinom.md)
  and
  [`pztnbinom2()`](https://janolefi.github.io/RTMBdist/reference/ztnbinom2.md),
  which returned `NaN` instead of 0 for quantiles below their support.

## RTMBdist 1.1.0

CRAN release: 2026-09-06

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
