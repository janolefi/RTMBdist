# Tests for the Sichel distribution

test_that("sichel passes standard discrete checks", {
  check_discrete_dist(dsichel, psichel, xs_int = c(0, 1, 4, 10, 30),
                      sum_support = 0:2000, mu = 4, sigma = 0.5, nu = 1.3)
  check_discrete_dist(dsichel, psichel, xs_int = c(0, 2, 8, 25, 80),
                      sum_support = 0:20000, mu = 5, sigma = 3, nu = -2)
})

test_that("dsichel and psichel match gamlss.dist", {
  # reference values from gamlss.dist::dSICHEL / pSICHEL, baked in: gamlss.dist
  # is scheduled for archival, so the test must not depend on it
  x <- c(0, 1, 4, 15, 45, 120, 2, 7)
  mu <- c(4, 4, 0.5, 10, 2, 25, 3, 6)
  s <- c(0.5, 3, 0.1, 0.7, 20, 0.05, 1, 2)
  nu <- c(1.3, -2, 0.7, -0.5, -1.2, 4, 0.2, -3.5)
  logd <- c(-2.5438928660815, -1.65858287991499, -6.09044949736467, -3.7524480663837,
            -9.01165220942829, -29.8854889368519, -1.78727496471604, -2.62682256396652)
  cdf <- c(0.0785599802700066, 0.317094219792606, 0.999650438413505, 0.814665496053839,
           0.99596315276136, 0.99999999999977, 0.56803417631402, 0.740657561981258)
  expect_equal(dsichel(x, mu, s, nu, log = TRUE), logd, tolerance = 1e-10)
  expect_equal(psichel(x, mu, s, nu), cdf, tolerance = 1e-10)
})

test_that("mu is the mean and the variance matches the closed form", {
  k <- 0:200000
  for (p in list(c(4, 0.5, 1.3), c(10, 2, -2))) {
    mu <- p[1]; sigma <- p[2]; nu <- p[3]
    d <- dsichel(k, mu, sigma, nu)
    c <- besselK(1 / sigma, nu + 1) / besselK(1 / sigma, nu)
    m <- sum(k * d)
    expect_equal(m, mu, tolerance = 1e-8)
    expect_equal(sum((k - m)^2 * d),
                 mu + mu^2 * (2 * sigma * (nu + 1) / c + 1 / c^2 - 1), tolerance = 1e-6)
  }
})

test_that("dsichel stays finite and accurate for large counts", {
  # The three-term recursion of the Bessel function carries over to the mass function,
  # P(x) = (x + nu - 1) / x * 2 mu sigma / (c s^2) P(x - 1) + (mu / (c s))^2 / (x (x - 1)) P(x - 2)
  # with s^2 = 1 + 2 mu sigma / c, which also checks the switch to the large-order
  # expansion of the Bessel function, at x + nu = 50
  for (p in list(c(5, 2, 1.3), c(2, 20, -0.5), c(40, 0.1, -3))) {
    mu <- p[1]; sigma <- p[2]; nu <- p[3]
    x <- 0:400
    ld <- dsichel(x, mu, sigma, nu, log = TRUE)
    expect_true(all(is.finite(ld)))
    c <- besselK(1 / sigma, nu + 1) / besselK(1 / sigma, nu)
    s2 <- 1 + 2 * mu * sigma / c
    i <- 3:length(x)
    xi <- x[i]
    rec <- log((xi + nu - 1) / xi * 2 * mu * sigma / (c * s2) * exp(ld[i - 1]) +
                 mu^2 / (c^2 * s2 * xi * (xi - 1)) * exp(ld[i - 2]))
    ok <- ld[i] > -700 # where the recursion itself does not underflow
    expect_equal(ld[i][ok], rec[ok], tolerance = 1e-9)
  }
})

test_that("dsichel gradients match finite differences, in all three parameters", {
  x <- c(0, 3, 12, 49, 52, 150)
  par <- c(mu = 5, sigma = 2, nu = 1.3)
  f <- function(p) sum(dsichel(x, p[1], p[2], p[3], log = TRUE))
  F <- RTMB::MakeTape(f, par)
  h <- 1e-5
  fd <- sapply(1:3, function(j) {
    e <- replace(numeric(3), j, h)
    (f(par + e) - f(par - e)) / (2 * h)
  })
  expect_equal(as.vector(F$jacobian(par)), fd, tolerance = 1e-6)
  expect_true(all(is.finite(F$jacfun()$jacobian(par))))
  check_ad_gradient(dsichel, rsichel, mu = 4, sigma = 0.5, nu = 1.3)
})

test_that("psichel is AD-compatible and gives OSA residuals", {
  q <- c(0, 1, 4, 10, 30)
  check_ad_cdf(psichel, dsichel, q, mu = 4, sigma = 0.5, nu = 1.3, .dq = FALSE)
  check_osa_cdf(dsichel, psichel, q, mu = 4, sigma = 0.5, nu = 1.3, .discrete = TRUE)
})

test_that("rsichel draws follow the mass function", {
  set.seed(3)
  x <- rsichel(2e4, 4, 0.5, 1.3)
  e <- c(dsichel(0:14, 4, 0.5, 1.3), 1 - sum(dsichel(0:14, 4, 0.5, 1.3)))
  o <- as.vector(table(factor(pmin(x, 15), levels = 0:15)))
  expect_gt(suppressWarnings(stats::chisq.test(o, p = e)$p.value), 0.01)
  expect_length(rsichel(5, c(2, 4), 0.5, 1), 5)
})

test_that("sichel recovers its parameters", {
  set.seed(11)
  x <- rsichel(20000, 6, 0.8, 1.5)
  f <- function(p) -sum(dsichel(x, exp(p[1]), exp(p[2]), p[3], log = TRUE))
  F <- RTMB::MakeTape(f, c(log(6), log(0.8), 1.5))
  o <- nlminb(c(0, 0, 0), F, function(p) as.vector(F$jacobian(p)))
  expect_equal(o$convergence, 0)
  expect_equal(c(exp(o$par[1:2]), o$par[3]), c(6, 0.8, 1.5), tolerance = 0.2)
})

test_that("sichel recycles, is zero below the support and rejects bad parameters", {
  expect_equal(dsichel(0:3, c(1, 2), 0.5, 1),
               c(dsichel(0, 1, 0.5, 1), dsichel(1, 2, 0.5, 1),
                 dsichel(2, 1, 0.5, 1), dsichel(3, 2, 0.5, 1)))
  expect_equal(dsichel(c(-3, -1), 4, 0.5, 1), c(0, 0))
  expect_error(dsichel(1, 0, 0.5, 1), "mu")
  expect_error(psichel(1, 4, -1, 1), "sigma")
  expect_error(rsichel(1, 4, 0, 1), "sigma")
})
