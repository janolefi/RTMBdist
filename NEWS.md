# RTMBdist 1.2.0

- Added the half-t (`dhalft()`) and half-Cauchy (`dhalfcauchy()`) distributions, each with matching `p`, `q` and `r` functions. These are the standard weakly informative priors for the standard deviation of a hierarchical model, so they are aimed squarely at models fitted by the Laplace approximation. Both distribution functions are differentiable, and `phalft()` differentiates with respect to `df` as well as `sigma`, so the degrees of freedom can be estimated rather than fixed; one-step-ahead residuals via `method = "cdf"` are supported. The half-Cauchy is the half-t with `df = 1`, and the half-normal is `dfoldnorm()` with `mu = 0`, which the half-t approaches as `df` grows.

- Added the Yule-Simon (`dyules()`) and Waring (`dwaring()`) distributions, two classical long-tailed count laws, each with matching `p` and `r` functions. Both are beta-geometric special cases of the beta-negative binomial, with `size` fixed at one, and that restriction is what gives them a closed-form distribution function where the general beta-negative binomial has none; one-step-ahead residuals via `method = "cdf"` are therefore available for these two. `dyules()` follows `VGAM` and is supported on the positive integers, while `dwaring()` follows the `WARING` family of `gamlss.dist`, starts at zero and has `mu` as its mean. The `YULE` family of `gamlss.dist` is `dwaring(x, mu, mu)`.

- Added the beta-negative binomial distribution (`dbnbinom()`) and its mean parameterisation (`dbnbinom2()`), each with a matching `r` function. It is to the negative binomial what the beta-binomial is to the binomial, and its extra beta layer gives a considerably heavier tail. In `dbnbinom2()`, which follows the `BNB` family of `gamlss.dist`, `mu` is exactly the mean; this is the more stable parameterisation to estimate in, because the original one has a long likelihood ridge along which `size` and `shape2` trade off. Like the beta-binomial, neither has a distribution function, so one-step-ahead residuals are not available.

- Added the extreme value distributions: the generalised extreme value distribution (`dgev()`), the generalised Pareto distribution (`dgpd()`) and the Frechet distribution (`dfrechet()`), each with matching `p`, `q` and `r` functions. The first two cover their three shape regimes with a single expression rather than a branch on the sign of `xi`, so the derivative with respect to the shape is exact at `xi = 0`, which is the usual starting value when the shape is estimated. Densities and distribution functions are differentiable, so simulation and one-step-ahead residuals via `method = "cdf"` are supported. Unlike `VGAM`, `evd` and `extraDistr`, `dgpd()` returns `1 / sigma` rather than zero at the threshold itself, matching `stats::dexp()` at zero.

- `RTMBdist` no longer masks anything in `stats`. The AD-compatible replacements for `stats::pt()`, `stats::plnorm()`, `stats::dgeom()` and `stats::pgeom()` are exported as `pt.ad()`, `plnorm.ad()`, `dgeom.ad()` and `pgeom.ad()`, and are reached through internal S4 generics that dispatch on the argument classes: plain numeric input goes to the `stats` versions, AD variables to the `.ad` versions. Previously `pt()` was exported with a reduced argument list, so `pt(q, df, lower.tail = FALSE)` failed for anyone who had loaded the package. As a side effect `lower.tail` and `log.p` now also work for `pt()` under automatic differentiation.

- Added the zero-inflated (`dzigeom()`), zero-truncated (`dztgeom()`) and hurdle (`dhgeom()`) geometric distributions. These build on `RTMB`'s AD-compatible negative binomial with `size = 1`, so `stats`' own `dgeom()` and friends are left untouched.

- Added the zero-inflated, zero-truncated and hurdle beta-binomial distributions (`dzibetabinom()`, `dztbetabinom()`, `dhbetabinom()`). Like `dbetabinom()` itself these have no distribution function, since the beta-binomial cdf has no closed form.

- The documentation of the continuous zero-inflated distributions now explains that, because the continuous part places no mass at zero, `zeroprob` is exactly the probability of a zero and zero-inflation coincides with a hurdle model; GAMLSS calls these zero-adjusted rather than zero-inflated.

- Added the hurdle (zero-altered) count distributions: Poisson (`dhpois()`), binomial (`dhbinom()`), negative binomial (`dhnbinom()`) and its mean parameterisation (`dhnbinom2()`), each with matching `p` and `r` functions. In a hurdle distribution the probability of a zero is a free parameter and the positive counts follow the corresponding zero-truncated distribution. Unlike zero-inflation, which can only add zeros, a hurdle model allows `zeroprob` to be smaller than the Poisson would give on its own. The density and distribution function are both differentiable, so simulation and one-step-ahead residuals via `method = "cdf"` are supported.

- Fixed `pztnbinom()` and `pztnbinom2()`, which returned `NaN` instead of 0 for quantiles below their support.

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
