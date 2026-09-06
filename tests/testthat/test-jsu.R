# Tests for the Johnson SU distribution, in both parameterisations

test_that("jsu passes standard distribution checks (mu=0, sigma=1, nu=0, tau=1)", {
  check_continuous_dist(
    dfun  = djsu,
    pfun  = pjsu,
    qfun  = qjsu,
    xs    = c(-2, -0.7, 0, 0.9, 3),
    mu = 0, sigma = 1, nu = 0, tau = 1
  )
})

test_that("jsu passes standard distribution checks (mu=1, sigma=2, nu=-1, tau=3)", {
  check_continuous_dist(
    dfun  = djsu,
    pfun  = pjsu,
    qfun  = qjsu,
    xs    = c(-1, 0.5, 1.5, 3, 6),
    mu = 1, sigma = 2, nu = -1, tau = 3
  )
})

test_that("jsu2 passes standard distribution checks (mu=3, sigma=2, nu=1, tau=3)", {
  check_continuous_dist(
    dfun  = djsu2,
    pfun  = pjsu2,
    qfun  = qjsu2,
    xs    = c(-1, 1, 3, 5, 8),
    mu = 3, sigma = 2, nu = 1, tau = 3
  )
})

test_that("jsu AD gradients have no NaN", {
  check_ad_gradient(djsu,  rjsu,  mu = 1, sigma = 2, nu = -1, tau = 3)
  check_ad_gradient(djsu2, rjsu2, mu = 1, sigma = 2, nu = -1, tau = 3)
})

test_that("jsu2 mu and sigma really are the mean and standard deviation", {
  for (pars in list(c(3, 2, 1, 3), c(0, 1, -2, 1.5), c(-5, 0.5, 0.5, 4))) {
    mu <- pars[1]; sigma <- pars[2]; nu <- pars[3]; tau <- pars[4]
    m1 <- integrate(function(z) z * djsu2(z, mu, sigma, nu, tau), -Inf, Inf)$value
    m2 <- integrate(function(z) z^2 * djsu2(z, mu, sigma, nu, tau), -Inf, Inf)$value
    expect_equal(m1, mu, tolerance = 1e-6)
    expect_equal(sqrt(m2 - m1^2), sigma, tolerance = 1e-6)
  }
})

test_that("jsu2 is the reparameterised jsu", {
  # djsu2(x, mu, sigma, nu, tau) == djsu(x, location, scale, -nu, tau)
  x <- c(-2, 0, 1.5, 4)
  mu <- 1; sigma <- 2; nu <- -1.5; tau <- 3
  rtau <- 1 / tau
  w <- exp(rtau^2); omega <- -nu * rtau
  cc <- 1 / sqrt(0.5 * (w - 1) * (w * cosh(2 * omega) + 1))
  expect_equal(
    djsu2(x, mu, sigma, nu, tau),
    djsu(x, mu + cc * sigma * sqrt(w) * sinh(omega), cc * sigma, -nu, tau)
  )
})

test_that("jsu tends to the normal distribution as tau grows", {
  x <- c(-2, -0.5, 0, 1, 2.5)
  # the moment parameterisation makes the limit exactly N(mu, sigma)
  expect_equal(djsu2(x, 0, 1, 0, 1e6), dnorm(x), tolerance = 1e-6)
  expect_equal(djsu2(x, 2, 3, 0, 1e8), dnorm(x, 2, 3), tolerance = 1e-6)
  expect_true(all(is.finite(djsu2(x, 0, 1, 0, 1e9))))
})

test_that("jsu quantile functions honour lower.tail and log.p", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  expect_equal(qjsu(p, 1, 2, -1, 3),  qjsu(1 - p, 1, 2, -1, 3, lower.tail = FALSE))
  expect_equal(qjsu(p, 1, 2, -1, 3),  qjsu(log(p), 1, 2, -1, 3, log.p = TRUE))
  expect_equal(qjsu2(p, 1, 2, -1, 3), qjsu2(1 - p, 1, 2, -1, 3, lower.tail = FALSE))
  expect_equal(qjsu2(p, 1, 2, -1, 3), qjsu2(log(p), 1, 2, -1, 3, log.p = TRUE))
})
