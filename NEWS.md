# RTMBdist 1.1.0

- Added the Johnson SU distribution in both the original parameterisation (`djsu()`, `pjsu()`, `qjsu()`, `rjsu()`) and the moment parameterisation (`djsu2()`, `pjsu2()`, `qjsu2()`, `rjsu2()`), a four-parameter distribution on the real line covering a wide range of skewness and kurtosis. In `djsu2()` the location and scale arguments are exactly the mean and standard deviation. The density and distribution function are both differentiable, so simulation and one-step-ahead residuals are supported. Unlike `gamlss.dist`, the reparameterisation stays finite for very large `tau`, where the distribution approaches the normal.

- Added references to the primary source for each distribution derived from `gamlss.dist`, and cross-links between related families. `pgenpois()` and friends previously had no references at all.

- Removed the dependency on `gamlss.dist`, which is scheduled for archival on CRAN. The quantile and random generation functions of the Box-Cox Cole-Green (`qbccg()`, `rbccg()`), Box-Cox *t* (`qbct()`, `rbct()`), Box-Cox power exponential (`qbcpe()`, `rbcpe()`), power exponential (`qpowerexp()`, `rpowerexp()`, `qpowerexp2()`, `rpowerexp2()`), Pareto (`qpareto()`, `rpareto()`), generalised Poisson (`pgenpois()`, `qgenpois()`, `rgenpois()`) and exponentially modified Gaussian (`qexgauss()`) distributions are now implemented natively. Results are unchanged except for the fixes below.

- Fixed argument recycling in `qbccg()`, `qbct()`, `qbcpe()`, `qgenpois()` and `pgenpois()`. Evaluating a single quantile against vectorised parameters previously collapsed the result to length one, and a parameter vector shorter than `x` silently truncated it. These functions now return one value per recycled argument tuple.

- Fixed `qbct()` and `rbct()`, which failed with an error when `mu`, `sigma` or `tau` was supplied as a vector.

- `qpareto()` now honours `lower.tail`, which was previously accepted but ignored.

- `log.p = TRUE` now works in `qpareto()`, `qbct()` and `qgenpois()`. These previously validated `p` before transforming it back from the log scale, so the argument could not be used.

- `qexgauss()` is substantially more accurate. It inverts `pexgauss()` with a tighter convergence tolerance, reducing the round-trip error `|p(q(p)) - p|` from roughly 1e-6 to roughly 1e-14.

- `pgenpois()` no longer falls back to the Poisson distribution for `phi < 1e-4` and evaluates the generalised Poisson distribution function across the whole parameter range.

- `pbetaprime()`, `pinvchisq()` and `pinvgamma()` are now differentiable, so one-step-ahead residuals via `method = "cdf"` are available for the beta prime, inverse chi-squared and inverse gamma distributions.

- Added the Bell distribution (`dbell()`, `pbell()`, `qbell()`, `rbell()`) and its mean parameterisation (`dbell2()`, `pbell2()`, `qbell2()`, `rbell2()`), a one-parameter distribution for overdispersed counts. Log Bell numbers are used throughout and cached, so the density stays finite well past the point at which the Bell numbers themselves overflow double precision. One-step-ahead residuals are supported via `method = "cdf"`.

- Added `lambertW()`, an AD-compatible implementation of the principal branch of the Lambert W function. Derivatives of all orders are exact, so it can be used in models fitted by Laplace approximation.

- Added the Conway-Maxwell-binomial distribution (`dcombinom()`, `pcombinom()`, `qcombinom()`, `rcombinom()`), a binomial generalisation with an additional dispersion parameter `nu`. The density and distribution function are differentiable with respect to `x`, so one-step-ahead residuals are supported.

# RTMBdist 0.1.0

- First release on CRAN.
