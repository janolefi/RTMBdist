# Tests for the half-t distribution

test_that("halft passes standard distribution checks", {
  check_continuous_dist(
    dfun = dhalft, pfun = phalft, qfun = qhalft,
    xs = c(0.1, 0.6, 1.5, 4, 20), lower = 0, upper = Inf,
    df = 4, sigma = 2
  )
  check_continuous_dist(
    dfun = dhalft, pfun = phalft, qfun = qhalft,
    xs = c(0.05, 0.4, 2, 15, 500), lower = 0, upper = Inf,
    df = 0.5, sigma = 1 # very heavy tail
  )
  check_continuous_dist(
    dfun = dhalft, pfun = phalft, qfun = qhalft,
    xs = c(0.2, 0.8, 1.6, 3, 6), lower = 0, upper = Inf,
    df = 30, sigma = 1.5 # close to half-normal
  )
})

test_that("halft AD gradient has no NaN", {
  check_ad_gradient(dhalft, rhalft, df = 4, sigma = 2)
  check_ad_gradient(dhalft, rhalft, df = 0.8, sigma = 1)
  check_ad_gradient(dhalft, rhalft, df = 50, sigma = 3)
})

test_that("dhalft is twice the scaled t density on the positive half line", {
  xs <- c(0, 0.3, 1, 4, 30)
  for (p in list(c(0.5, 1), c(3, 2), c(20, 0.4))) {
    expect_equal(dhalft(xs, p[1], p[2]), 2 * stats::dt(xs / p[2], p[1]) / p[2])
    expect_equal(phalft(xs, p[1], p[2]), 2 * stats::pt(xs / p[2], p[1]) - 1)
  }
})

test_that("halft mass sits on the closed positive half line", {
  expect_equal(dhalft(c(-4, -1, -1e-9), 3, 2), c(0, 0, 0))
  expect_equal(dhalft(c(-4, -1), 3, 2, log = TRUE), c(-Inf, -Inf))
  expect_equal(phalft(c(-4, -1), 3, 2), c(0, 0))
  # the origin belongs to the support and carries a positive density
  expect_equal(dhalft(0, 3, 2), 2 * stats::dt(0, 3) / 2)
  expect_equal(phalft(0, 3, 2), 0)
  expect_equal(qhalft(0, 3, 2), 0)
  expect_equal(qhalft(1, 3, 2), Inf)
})

test_that("halft has the half-Cauchy and half-normal as special cases", {
  xs <- c(0, 0.2, 1, 3, 10)
  expect_equal(dhalft(xs, 1, 2), dhalfcauchy(xs, 2))
  expect_equal(phalft(xs, 1, 2), phalfcauchy(xs, 2))
  # as df grows the half-t approaches the folded normal with mu = 0
  expect_equal(dhalft(xs, 1e7, 2), dfoldnorm(xs, 0, 2), tolerance = 1e-6)
})

test_that("halft moments match the closed forms", {
  for (d in c(3, 5, 10)) {
    m <- 2 * 2 * sqrt(d) * gamma((d + 1) / 2) / (sqrt(pi) * (d - 1) * gamma(d / 2))
    expect_equal(integrate(function(z) z * dhalft(z, d, 2), 0, Inf)$value, m,
                 tolerance = 1e-6)
    expect_equal(integrate(function(z) (z - m)^2 * dhalft(z, d, 2), 0, Inf)$value,
                 4 * d / (d - 2) - m^2, tolerance = 1e-5)
  }
})

test_that("phalft is differentiable in df as well as sigma", {
  # pt() dispatches to the AD implementation, so the degrees of freedom can be
  # estimated rather than fixed
  F <- RTMB::MakeTape(function(p) sum(phalft(c(0.5, 3), p[1], p[2])), c(4, 2))
  g <- as.vector(F$jacobian(c(4, 2)))
  num <- sapply(1:2, function(i) {
    h <- 1e-6; u <- l <- c(4, 2); u[i] <- u[i] + h; l[i] <- l[i] - h
    (sum(phalft(c(0.5, 3), u[1], u[2])) - sum(phalft(c(0.5, 3), l[1], l[2]))) / (2 * h)
  })
  expect_false(any(is.nan(g)))
  expect_equal(g, num, tolerance = 1e-6)
})

test_that("qhalft honours lower.tail and log.p", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  expect_equal(qhalft(p, 3, 2), qhalft(1 - p, 3, 2, lower.tail = FALSE))
  expect_equal(qhalft(p, 3, 2), qhalft(log(p), 3, 2, log.p = TRUE))
})

test_that("dhalft rejects non-positive parameters", {
  expect_error(dhalft(1, 0, 1), "df")
  expect_error(dhalft(1, 3, -1), "sigma")
  expect_error(phalft(1, -2, 1), "df")
})
