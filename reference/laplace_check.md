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

  Number of decimal places for the diagnostic values.

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
# Chicken weight example; taken from RTMB Introduction vignette
data(ChickWeight)

parameters <- list(
  mua=0,          ## Mean slope
  sda=1,          ## Std of slopes
  mub=0,          ## Mean intercept
  sdb=1,          ## Std of intercepts
  sdeps=1,        ## Residual Std
  a=rep(0, 50),   ## Random slope by chick
  b=rep(0, 50)    ## Random intercept by chick
)

jnll <- function(parms) {
  getAll(ChickWeight, parms, warn=FALSE)
  ## Optional (enables extra RTMB features)
  weight <- OBS(weight)
  ## Initialize joint negative log likelihood
  nll <- 0
  ## Random slopes
  nll <- nll - sum(dnorm(a, mean=mua, sd=sda, log=TRUE))
  ## Random intercepts
  nll <- nll - sum(dnorm(b, mean=mub, sd=sdb, log=TRUE))
  ## Data
  predWeight <- a[Chick] * Time + b[Chick]
  nll <- nll - sum(dnorm(weight, predWeight, sd=sdeps, log=TRUE))
  ## Get predicted weight uncertainties
  ADREPORT(predWeight)
  ## Return
  nll
}

obj <- MakeADFun(jnll, parameters, random=c("a", "b"), silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)

chk <- laplace_check(obj)
chk
#> Laplace approximation check
#> ===========================
#> Monte Carlo samples:         1000
#> 
#> Log marginal likelihood
#>   Laplace:                   -2446.99
#>   Importance sampling:       -2446.99   (unbiased; should match Laplace)
#>   Bias (IS - Laplace):       0.0000   (log scale; 0 = exact)
#>     likelihood ratio:        1.000   (exp(bias); 1 = no error)
#>     relative to log-lik:     1.9e-16   (0 = exact)
#> 
#> Diagnostics
#>   SD of log-weights:         0.0000   (0 = exact; grows with non-Gaussianity and dimension)
#>   Effective sample size:     100.0%   (100% = ideal; low may reflect high dimension)
#> 
#> The Laplace approximation appears accurate.
# Laplace approximation exact here: linear Gaussian-Gaussian example
```
