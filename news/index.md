# Changelog

## RTMBdist (development version)

- Ongoing development version.

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
