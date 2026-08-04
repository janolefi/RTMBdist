# Assess the accuracy of the Laplace approximation

Diagnoses how well the Laplace approximation used by `TMB`/`RTMB`
approximates the marginal likelihood, via the importance-sampling check
of evaluating the integrand relative to the fitted Gaussian.

## Usage

``` r
laplace_check(obj, nSamples = 1000)

# S3 method for class 'laplace_check'
print(x, digits = 4, ...)
```

## Arguments

- obj:

  A fitted `TMB`/`RTMB` object (as returned by `MakeADFun`) with random
  effects. The model is assumed to have been optimised, so that
  `obj$env$last.par.best` holds the maximum-likelihood parameters and
  the corresponding mode of the latent variables.

- nSamples:

  Number of Monte Carlo samples drawn from the Gaussian approximation.
  Default `1000`.

- x:

  An object of class `"laplace_check"`.

- digits:

  Number of significant digits for printed values.

- ...:

  Unused.

## Value

An object of class `"laplace_check"`: a list with the log-weights
`logw`, their standard deviation `sd_logw`, the Laplace and
importance-sampling estimates of the log marginal likelihood
(`lZ_laplace`, `lZ_is`), their difference `log_bias`, and the relative
effective sample size `ess_ratio`.

## Details

For fixed parameters \\\hat\theta\\, the marginal likelihood is \\Z =
\int f(y\mid x) f(x)\\ dx\\. Writing the importance weight
\$\$\tilde\ell(x) = f(y\mid x) f(x) / N(x\mid \hat x, H^{-1}),\$\$ the
Laplace approximation equals \\\tilde\ell(\hat x)\\, while \\E\_{x\sim
N(\hat x, H^{-1})}\[\tilde\ell(x)\] = Z\\ exactly. The two coincide if
and only if the joint negative log-likelihood is quadratic in the latent
variables. The spread of the log-weights therefore measures the
departure from this ideal: a standard deviation near zero indicates a
near-exact approximation.

## Examples

``` r
# currently no example
```
