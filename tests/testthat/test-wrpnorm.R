# Tests for the wrapped normal distribution (circular). The distribution function
# has its origin at mu - pi, so it is checked on [mu - pi, mu + pi).

test_that("wrpnorm passes standard distribution checks on both sides of sigma = 2", {
  check_continuous_dist(dwrpnorm, pwrpnorm, NULL, xs = c(-3, -1, 0, 0.5, 3),
                        lower = -pi, upper = pi, mu = 0, sigma = 0.5)
  check_continuous_dist(dwrpnorm, pwrpnorm, NULL, xs = c(-2, 0, 1, 1.2, 3.5),
                        lower = 1 - pi, upper = 1 + pi, mu = 1, sigma = 3)
})

test_that("dwrpnorm matches the sum over many windings", {
  # circular::dwrappednormal is off by 1e-10 for sigma = 5, so it is not used here
  x <- seq(-pi, pi, length.out = 13)
  for (sigma in c(0.1, 0.7, 1.9, 2, 2.1, 5, 20)) {
    d <- sapply(x, function(xi) sum(dnorm(xi - 1 + 2 * pi * (-300:300), 0, sigma)))
    expect_equal(dwrpnorm(x, 1, sigma), d, tolerance = 1e-12)
  }
})

test_that("the two series agree where they switch at sigma = 2", {
  x <- seq(-pi, pi, length.out = 9)
  expect_equal(dwrpnorm(x, 0.5, 2 - 1e-12), dwrpnorm(x, 0.5, 2), tolerance = 1e-12)
  expect_equal(pwrpnorm(x, 0.5, 2 - 1e-12), pwrpnorm(x, 0.5, 2), tolerance = 1e-12)
})

test_that("wrpnorm is the normal for small sigma and the uniform for large sigma", {
  x <- c(-0.3, 0, 0.2)
  expect_equal(dwrpnorm(x, 0, 0.1), dnorm(x, 0, 0.1))
  expect_equal(dwrpnorm(x, 0, 50), rep(1 / (2 * pi), 3))
  expect_equal(pwrpnorm(x, 0, 50), (x + pi) / (2 * pi))
})

test_that("dwrpnorm and pwrpnorm are periodic, and the origin is at mu - pi", {
  xs <- c(-3, -1, 0.5, 2.5)
  for (sigma in c(0.6, 3)) {
    expect_equal(dwrpnorm(xs + 4 * pi, 0.7, sigma), dwrpnorm(xs, 0.7, sigma))
    expect_equal(pwrpnorm(xs - 2 * pi, 0.7, sigma), pwrpnorm(xs, 0.7, sigma))
    expect_equal(pwrpnorm(0.7, 0.7, sigma), 0.5)
  }
})

test_that("pwrpnorm with from integrates the density from that origin", {
  for (sigma in c(0.6, 3)) {
    for (from in c(-pi, 0.4, 10)) {
      xs <- from + c(0.3, 2.5, 6)
      num <- sapply(xs, function(x)
        integrate(dwrpnorm, from, x, mu = 2, sigma = sigma, rel.tol = 1e-10)$value)
      expect_equal(pwrpnorm(xs, 2, sigma, from = from), num, tolerance = 1e-8)
    }
  }
  expect_equal(pwrpnorm(1, 1, 0.5, from = 1 - pi, lower.tail = FALSE),
               pwrpnorm(1, 1, 0.5, lower.tail = FALSE))
})

test_that("the upper tail is accurate far from the mean", {
  # 1 - p would be 0 here
  expect_equal(pwrpnorm(0.5, 0, 0.05, lower.tail = FALSE), pnorm(-0.5 / 0.05))
})

test_that("wrpnorm AD derivatives are correct", {
  check_ad_gradient(dwrpnorm, rwrpnorm, mu = 0.3, sigma = 0.8)
  check_ad_gradient(dwrpnorm, rwrpnorm, mu = 0.3, sigma = 2.5)
  xs <- c(-3, -1, 0, 1, 3)
  check_ad_cdf(pwrpnorm, dwrpnorm, xs, mu = 0.3, sigma = 0.8)
  check_ad_cdf(pwrpnorm, dwrpnorm, xs, mu = 0.3, sigma = 2.5)
  check_ad_cdf_hessian(pwrpnorm, xs, mu = 0.3, sigma = 0.8)
})

test_that("pwrpnorm with a fixed origin stays correct when the tape is re-evaluated", {
  q <- c(-2, 0, 2)
  F <- RTMB::MakeTape(function(m) pwrpnorm(q, m, 0.5, from = -pi), 0)
  expect_equal(F(2.5), pwrpnorm(q, 2.5, 0.5, from = -pi), tolerance = 1e-10)
})

test_that("wrpnorm supports OSA residuals, with the circle cut at -pi", {
  set.seed(1)
  check_osa_cdf(dwrpnorm, function(q, mu, sigma) pwrpnorm(q, mu, sigma, from = -pi),
                rwrpnorm(30, 1, 0.8), mu = 1, sigma = 0.8)
})

test_that("rwrpnorm draws from the wrapped normal", {
  set.seed(1)
  x <- rwrpnorm(20000, mu = 2.5, sigma = 0.9)
  expect_true(all(x >= -pi & x < pi))
  expect_gt(ks.test(pwrpnorm(x, 2.5, 0.9, from = -pi), "punif")$p.value, 0.01)
  expect_length(rwrpnorm(4, mu = 1:4, sigma = 0.5), 4)
})

test_that("wrpnorm recycles and rejects non-positive sigma", {
  expect_equal(dwrpnorm(0.5, 0, c(0.5, 3)), c(dwrpnorm(0.5, 0, 0.5), dwrpnorm(0.5, 0, 3)))
  expect_equal(pwrpnorm(0.5, 0, c(0.5, 3)), c(pwrpnorm(0.5, 0, 0.5), pwrpnorm(0.5, 0, 3)))
  expect_error(dwrpnorm(1, 0, 0), "sigma")
  expect_error(pwrpnorm(1, 0, -1), "sigma")
  expect_error(rwrpnorm(1, 0, 0), "sigma")
})
