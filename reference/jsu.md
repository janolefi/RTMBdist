# Johnson SU distribution (JSU)

Density, distribution function, quantile function, and random generation
for the Johnson SU distribution, in the original and in the moment
parameterisation.

## Usage

``` r
djsu(x, mu = 0, sigma = 1, nu = 0, tau = 1, log = FALSE)

pjsu(q, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE)

qjsu(p, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE)

rjsu(n, mu = 0, sigma = 1, nu = 0, tau = 1)

djsu2(x, mu = 0, sigma = 1, nu = 0, tau = 1, log = FALSE)

pjsu2(q, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE)

qjsu2(p, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE)

rjsu2(n, mu = 0, sigma = 1, nu = 0, tau = 1)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter for `djsu`; the mean for `djsu2`.

- sigma:

  scale parameter for `djsu`; the standard deviation for `djsu2`. Must
  be positive.

- nu:

  skewness parameter (real). Positive \\\nu\\ gives left skewness in
  `djsu` and right skewness in `djsu2`.

- tau:

  kurtosis parameter, must be positive. Large \\\tau\\ approaches the
  normal distribution.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`djsu` gives the density, `pjsu` gives the distribution function, `qjsu`
gives the quantile function, and `rjsu` generates random deviates.
`djsu2`, `pjsu2`, `qjsu2` and `rjsu2` are the corresponding functions
for the moment parameterisation.

## Details

The Johnson SU distribution is a four-parameter continuous distribution
on the whole real line, obtained by applying the transformation \$\$Z =
\nu + \tau \\\mathrm{asinh}\\\left(\frac{x-\mu}{\sigma}\right) \sim
N(0,1)\$\$ to a standard normal variable. It covers a wide range of
skewness and kurtosis combinations and is a common alternative to the
Box-Cox families for data that are not restricted to be positive.

`djsu` uses the original parameterisation, in which \\\mu\\ and
\\\sigma\\ are a location and a scale parameter, \\\nu\\ controls
skewness and \\\tau\\ controls kurtosis. The density is \$\$f(x; \mu,
\sigma, \nu, \tau) = \frac{\tau}{\sigma \sqrt{2\pi}}
\frac{1}{\sqrt{z^2 + 1}} \exp\\\left(-\frac{r^2}{2}\right),\$\$ where
\\z = (x-\mu)/\sigma\\ and \\r = \nu + \tau\\\mathrm{asinh}(z)\\. Here
\\\mu\\ is *not* the mean and \\\sigma\\ is *not* the standard
deviation.

`djsu2` uses the moment parameterisation, in which \\\mu\\ **is** the
mean and \\\sigma\\ **is** the standard deviation of the distribution,
for any admissible \\\nu\\ and \\\tau\\. This is usually the more
convenient parameterisation for regression modelling, because the
location and scale parameters keep their interpretation as \\\nu\\ and
\\\tau\\ change. Writing \\\omega = -\nu/\tau\\ and \\w =
\exp(\tau^{-2})\\, it is obtained from the original parameterisation by
the reparameterisation \$\$c = \left\[\tfrac{1}{2}(w-1)\left(w
\cosh(2\omega) + 1\right)\right\]^{-1/2}, \qquad \sigma^\* = c\\\sigma,
\qquad \mu^\* = \mu + c\\\sigma\sqrt{w}\\\sinh(\omega),\$\$ so that
`djsu2(x, mu, sigma, nu, tau)` equals
`djsu(x, `\\\mu^\*\\`, `\\\sigma^\*\\`, -nu, tau)`.

These correspond to the `JSUo` and `JSU` families of the `gamlss.dist`
package respectively; see Chapter 18 of Rigby et al. (2019). Note that
the sign convention for \\\nu\\ differs between the two
parameterisations, exactly as it does in `gamlss.dist`.

All four `d` and `p` functions are compatible with automatic
differentiation by `RTMB`, so both simulation and one-step-ahead
residuals are supported.

## References

Johnson, N. L. (1954). Systems of frequency curves derived from the
first law of Laplace. Trabajos de Estadistica, 5, 283-291.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## See also

[bcpe](https://janolefi.github.io/RTMBdist/reference/bcpe.md),
[bct](https://janolefi.github.io/RTMBdist/reference/bct.md),
[skewt](https://janolefi.github.io/RTMBdist/reference/skewt.md),
[skewnorm](https://janolefi.github.io/RTMBdist/reference/skewnorm.md)

## Examples

``` r
set.seed(123)
# original parameterisation
x <- rjsu(5, mu = 0, sigma = 1, nu = -1, tau = 2)
d <- djsu(x, mu = 0, sigma = 1, nu = -1, tau = 2)
p <- pjsu(x, mu = 0, sigma = 1, nu = -1, tau = 2)
q <- qjsu(p, mu = 0, sigma = 1, nu = -1, tau = 2)

# moment parameterisation: mu is the mean, sigma the standard deviation
y <- rjsu2(1000, mu = 3, sigma = 2, nu = 1, tau = 3)
c(mean = mean(y), sd = sd(y))
#>     mean       sd 
#> 2.972965 1.979446 
```
