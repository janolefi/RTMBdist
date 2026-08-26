# RTMBdist (development version)

- Ongoing development version.

- Added the Bell distribution (`dbell()`, `pbell()`, `qbell()`, `rbell()`) and its mean parameterisation (`dbell2()`, `pbell2()`, `qbell2()`, `rbell2()`), a one-parameter distribution for overdispersed counts. Log Bell numbers are used throughout and cached, so the density stays finite well past the point at which the Bell numbers themselves overflow double precision. One-step-ahead residuals are supported via `method = "cdf"`.

- Added `lambertW()`, an AD-compatible implementation of the principal branch of the Lambert W function. Derivatives of all orders are exact, so it can be used in models fitted by Laplace approximation.

- Added the Conway-Maxwell-binomial distribution (`dcombinom()`, `pcombinom()`, `qcombinom()`, `rcombinom()`), a binomial generalisation with an additional dispersion parameter `nu`. The density and distribution function are differentiable with respect to `x`, so one-step-ahead residuals are supported.

# RTMBdist 0.1.0

- First release on CRAN.
