# List of distributions

## Continuous distributions

- [`bccg(mu, sigma, nu)`](https://janoleko.github.io/RTMBdist/reference/bccg.md):
  Box-Cox Cole and Green distribution parameterised by location `mu`,
  scale `sigma`, and skewness `nu`

- [`bcpe(mu, sigma, nu, tau)`](https://janoleko.github.io/RTMBdist/reference/bcpe.md):
  Box-Cox power exponential distribution parameterised by location `mu`,
  scale `sigma`, `nu`, and `tau`

- [`bct(mu, sigma, nu, tau)`](https://janoleko.github.io/RTMBdist/reference/bct.md):
  Box-Cox t-distribution parameterised by location `mu`, scale `sigma`,
  skewness `nu`, and degrees of freedom `tau`

- [`beta2(mu, phi)`](https://janoleko.github.io/RTMBdist/reference/beta2.md):
  Beta distribution reparameterised by mean `mu` and precision `phi`

- [`exgauss(mu, sigma, lambda)`](https://janoleko.github.io/RTMBdist/reference/exgauss.md):
  Exponentially modified Gaussian distribution parameterised by location
  `mu`, scale `sigma`, and rate `lambda`

- [`foldnorm(mu, sigma)`](https://janoleko.github.io/RTMBdist/reference/foldnorm.md):
  Folded normal distribution parameterised by location `mu` and scale
  `sigma`

- [`gamma2(mean, sd)`](https://janoleko.github.io/RTMBdist/reference/gamma2.md):
  Gamma distribution reparameterised by mean and standard deviation

- [`gengamma(mu, sigma, nu)`](https://janoleko.github.io/RTMBdist/reference/gengamma.md):
  Generalised gamma distribution parameterised by location `mu`, scale
  `sigma`, and skewness `nu`

- [`gumbel(location, scale)`](https://janoleko.github.io/RTMBdist/reference/gumbel.md):
  Gumbel distribution parameterised by location and scale

- [`invgauss(mean, shape)`](https://janoleko.github.io/RTMBdist/reference/invgauss.md):
  Inverse Gaussian distribution parameterised by mean and shape

- [`laplace(mu, b)`](https://janoleko.github.io/RTMBdist/reference/laplace.md):
  Laplace distribution parameterised by location `mu` and scale `b`

- [`oibeta(shape1, shape2, oneprob)`](https://janoleko.github.io/RTMBdist/reference/oibeta.md):
  One-inflated beta distribution parameterised by shape parameters
  `shape1`, `shape2` and one-probability `oneprob`

- [`oibeta2(mu, phi, oneprob)`](https://janoleko.github.io/RTMBdist/reference/oibeta2.md):
  One-inflated beta distribution reparameterised by mean `mu`, precision
  `phi`, and one-probability `oneprob`

- [`pareto(mu)`](https://janoleko.github.io/RTMBdist/reference/pareto.md):
  Pareto distribution parameterised by `mu`

- [`powerexp(mu, sigma, nu)`](https://janoleko.github.io/RTMBdist/reference/powerexp.md):
  Power exponential distribution parameterised by mean `mu`, standard
  deviation `sigma` and shape `nu`

- [`powerexp2(mu, sigma, nu)`](https://janoleko.github.io/RTMBdist/reference/powerexp.md):
  Power exponential distribution reparameterised by location `mu`, scale
  `sigma` and shape `nu`

- [`pgweibull(scale, shape, powershape)`](https://janoleko.github.io/RTMBdist/reference/pgweibull.md):
  Power generalised Weibull distribution parameterised by `scale`,
  `shape` and `powershape`

- [`skewnorm(xi, omega, alpha)`](https://janoleko.github.io/RTMBdist/reference/skewnorm.md):
  Skew normal distribution parameterised by location `xi`, scale `omega`
  and skewness `alpha`

- [`skewnorm2(mean, sd, alpha)`](https://janoleko.github.io/RTMBdist/reference/skewnorm2.md):
  Skew normal distribution reparameterised by mean, standard deviation
  and skewness `alpha`

- [`skewt(mu, sigma, skew, df)`](https://janoleko.github.io/RTMBdist/reference/skewt.md):
  Skew t-distribution parameterised by location `mu`, scale `sigma`,
  skewness `skew` and degrees of freedom `df`

- [`skewt2(mean, sd, skew, df)`](https://janoleko.github.io/RTMBdist/reference/skewt2.md):
  Skew t-distribution reparameterised by mean, standard deviation,
  skewness `skew` and degrees of freedom `df`

- [`truncnorm(mean, sd, min, max)`](https://janoleko.github.io/RTMBdist/reference/truncnorm.md):
  Truncated normal distribution parameterised by mean, standard
  deviation, lower bound `min` and upper bound `max`

- [`trunct(df, min, max)`](https://janoleko.github.io/RTMBdist/reference/trunct.md):
  Truncated t-distribution parameterised by degrees of freedom `df`,
  lower bound `min` and upper bound `max`

- [`trunct2(df, mu, sigma, min, max)`](https://janoleko.github.io/RTMBdist/reference/trunct.md):
  Truncated t-distribution parameterised location `mu`, scale `sigma`,
  degrees of freedom `df`, lower bound `min` and upper bound `max`

- [`t2(mu, sigma, df)`](https://janoleko.github.io/RTMBdist/reference/t2.md):
  Non-central and scaled t-distribution parameterised by location `mu`,
  scale `sigma` and degrees of freedom `df`

- [`vm(mu, kappa)`](https://janoleko.github.io/RTMBdist/reference/vm.md):
  Von Mises distribution parameterised by mean direction `mu` and
  concentration `kappa`

- [`wrpcauchy(mu, rho)`](https://janoleko.github.io/RTMBdist/reference/wrpcauchy.md):
  Wrapped Cauchy distribution parameterised by mean direction `mu` and
  concentration `rho`

- [`zibeta(shape1, shape2, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/zibeta.md):
  Zero-inflated beta distribution parameterised by shape parameters
  `shape1`, `shape2` and zero-probability `zeroprob`

- [`zibeta2(mu, phi, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/zibeta2.md):
  Zero-inflated beta distribution reparameterised by mean `mu`,
  precision `phi`, and zero-probability `zeroprob`

- [`zigamma(shape, scale, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/zigamma.md):
  Zero-inflated gamma distribution parameterised by shape and scale,
  with a zero-probability `zeroprob`

- [`zigamma2(mean, sd, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/zigamma2.md):
  Zero-inflated gamma distribution reparameterised by mean, standard
  deviation and zero-probability `zeroprob`

- [`ziinvgauss(mean, shape, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/ziinvgauss.md):
  Zero-inflated inverse Gaussian distribution parameterised by mean,
  shape and zero-probability `zeroprob`

- [`zilnorm(meanlog, sdlog, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/zilnorm.md):
  Zero-inflated log normal distribution parameterised by meanlog, sdlog
  and zero-probability `zeroprob`

- [`zoibeta(shape1, shape2, zeroprob, oneprob)`](https://janoleko.github.io/RTMBdist/reference/zoibeta.md):
  Zero- and one-inflated beta distribution parameterised by shape
  parameters `shape1`, `shape2`, zero-probability `zeroprob` and
  one-probability `oneprob`

- [`zoibeta2(mu, phi, zeroprob, oneprob)`](https://janoleko.github.io/RTMBdist/reference/zoibeta2.md):
  Zero- and one-inflated beta distribution reparameterised by mean `mu`,
  precision `phi`, zero-probability `zeroprob` and one-probability
  `oneprob`

## Discrete distributions

- [`betabinom(size, shape1, shape2)`](https://janoleko.github.io/RTMBdist/reference/betabinom.md):
  Beta-binomial distribution parameterised by size `size`, shape
  parameters `shape1` and `shape2`

- [`genpois(lambda, phi)`](https://janoleko.github.io/RTMBdist/reference/genpois.md):
  Generalised Poisson distribution parameterised by mean `lambda` and
  dispersion `phi`

- [`nbinom2(mu, size)`](https://janoleko.github.io/RTMBdist/reference/nbinom2.md):
  Negative binomial distribution reparameterised by mean `mu` and size
  `size`

- [`zibinom(size, prob, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/zibinom.md):
  Zero-inflated binomial distribution parameterised by size `size`,
  success probability `prob` and zero-probability `zeroprob`

- [`zinbinom(size, prob, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/zinbinom.md):
  Zero-inflated negative binomial distribution parameterised by size
  `size`, success probability `prob` and zero-probability `zeroprob`

- [`zinbinom2(mu, size, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/zinbinom2.md):
  Zero-inflated negative binomial distribution reparameterised by mean
  `mu`, size `size` and zero-probability `zeroprob`

- [`zipois(lambda, zeroprob)`](https://janoleko.github.io/RTMBdist/reference/zipois.md):
  Zero-inflated Poisson distribution parameterised by rate `lambda` and
  zero-probability `zeroprob`

- [`ztbinom(size, prob)`](https://janoleko.github.io/RTMBdist/reference/ztbinom.md):
  Zero-truncated binomial distribution parameterised by size `size` and
  success probability `prob`

- [`ztnbinom(size, prob)`](https://janoleko.github.io/RTMBdist/reference/ztnbinom.md):
  Zero-truncated negative binomial distribution parameterised by size
  `size` and success probability `prob`

- [`ztnbinom2(mu, size)`](https://janoleko.github.io/RTMBdist/reference/ztnbinom2.md):
  Zero-truncated negative binomial distribution reparameterised by mean
  `mu` and size `size`

- [`ztpois(lambda)`](https://janoleko.github.io/RTMBdist/reference/ztpois.md):
  Zero-truncated Poisson distribution parameterised by rate `lambda`

## Multivariate distributions

- [`dirichlet(alpha)`](https://janoleko.github.io/RTMBdist/reference/dirichlet.md):
  Dirichlet distribution parameterised by concentration parameter vector
  `alpha`

- [`dirmult(size, alpha)`](https://janoleko.github.io/RTMBdist/reference/dirmult.md):
  Dirichlet-multinomial distribution parameterised by `size` and
  concentration parameters `alpha`

- [`mvt(mu, Sigma, df)`](https://janoleko.github.io/RTMBdist/reference/mvt.md):
  Multivariate t-distribution parameterised by location `mu`, scale
  matrix `Sigma` and degrees of freedom `df`

- [`vmf(mu, kappa)`](https://janoleko.github.io/RTMBdist/reference/vmf.md):
  Multivariate von Mises-Fisher distribution parameterised by unit mean
  vector `mu` and concentration `kappa`

- [`vmf2(theta)`](https://janoleko.github.io/RTMBdist/reference/vmf2.md):
  Multivariate von Mises-Fisher distribution parameterised by parameter
  `theta` equal to unit mean vector `mu` times concentration scalar
  `kappa`

- [`wishart(nu, Sigma)`](https://janoleko.github.io/RTMBdist/reference/wishart.md):
  Wishart distribution parameterised by degrees of freedom `nu` and
  scale matrix `Sigma`

## Copulas

Bivariate copulas can be implemented in a modular way using the
[`dcopula`](https://janoleko.github.io/RTMBdist/reference/dcopula.md)
function together with one of the copula constructors below. Available
copula constructors are:

- [`cgaussian(rho)`](https://janoleko.github.io/RTMBdist/reference/cgaussian.md)
  (Gaussian copula)
- [`cclayton(theta)`](https://janoleko.github.io/RTMBdist/reference/cclayton.md)
  (Clayton copula)
- [`cgumbel(theta)`](https://janoleko.github.io/RTMBdist/reference/cgumbel.md)
  (Gumbel copula)
- [`cfrank(theta)`](https://janoleko.github.io/RTMBdist/reference/cfrank.md)
  (Frank copula)

For bivariate copulas with *discrete* margins, use the
[`ddcopula`](https://janoleko.github.io/RTMBdist/reference/ddcopula.md)
function instead. In this case, instead of copula *densities*, copula
*CDFs* are needed. The available constructors for this are:

- [`Cclayton`](https://janoleko.github.io/RTMBdist/reference/Cclayton.md)
  (Clayton copula CDF)
- [`Cgumbel`](https://janoleko.github.io/RTMBdist/reference/Cgumbel.md)
  (Gumbel copula CDF)
- [`Cfrank`](https://janoleko.github.io/RTMBdist/reference/Cfrank.md)
  (Frank copula CDF)

Multivariate copulas are also possible using the
[`dmvcopula`](https://janoleko.github.io/RTMBdist/reference/dmvcopula.md)
function together with one of the multivariate copula constructors
below. Currently, only the multivariate Gaussian copula is implemented
in two ways:

- [`cmvgauss(R)`](https://janoleko.github.io/RTMBdist/reference/cmvgauss.md)
  (multivariate Gaussian copula parameterised by a correlation matrix)
- [`cgmrf(Q)`](https://janoleko.github.io/RTMBdist/reference/cgmrf.md)
  (multivariate Gaussian copula parameterised by an inverse correlation
  matrix)
