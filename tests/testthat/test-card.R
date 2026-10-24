# Tests for the cardioid distribution (circular). The distribution function has
# its origin at mu - pi, so it is checked on [mu - pi, mu + pi).

test_that("card passes standard distribution checks", {
  check_continuous_dist(dcard, pcard, NULL, xs = c(-3, -1, 0, 0.5, 3),
                        lower = -pi, upper = pi, mu = 0, rho = 0.3)
  check_continuous_dist(dcard, pcard, NULL, xs = c(-2, 0, 1, 1.2, 3.5),
                        lower = 1 - pi, upper = 1 + pi, mu = 1, rho = 0.5)
})

test_that("card is the circular uniform for rho = 0", {
  xs <- c(-3, -1, 0, 2, 3)
  expect_equal(dcard(xs, 1, 0), rep(1 / (2 * pi), 5))
  expect_equal(pcard(xs, 0, 0), (xs + pi) / (2 * pi))
})

test_that("pcard is periodic, has its origin at mu - pi, and honours from", {
  xs <- c(-3, -1, 0.5, 2.5)
  expect_equal(pcard(xs - 2 * pi, 0.7, 0.4), pcard(xs, 0.7, 0.4))
  expect_equal(pcard(0.7, 0.7, 0.4), 0.5)
  for (from in c(-pi, 0.4, 10)) {
    xs <- from + c(0.3, 2.5, 6)
    num <- sapply(xs, function(x) integrate(dcard, from, x, mu = 2, rho = 0.4)$value)
    expect_equal(pcard(xs, 2, 0.4, from = from), num, tolerance = 1e-8)
  }
})

test_that("card AD derivatives are correct", {
  check_ad_gradient(dcard, rcard, mu = 0.3, rho = 0.25)
  xs <- c(-3, -1, 0, 1, 3)
  check_ad_cdf(pcard, dcard, xs, mu = 0.3, rho = 0.25)
})

test_that("card supports OSA residuals, with the circle cut at -pi", {
  set.seed(1)
  check_osa_cdf(dcard, function(q, mu, rho) pcard(q, mu, rho, from = -pi),
                rcard(30, 1, 0.3), mu = 1, rho = 0.3)
})

test_that("rcard draws from the cardioid", {
  set.seed(1)
  x <- rcard(20000, mu = 2.5, rho = 0.4)
  expect_true(all(x >= -pi & x < pi))
  expect_gt(ks.test(pcard(x, 2.5, 0.4, from = -pi), "punif")$p.value, 0.01)
  # the mean resultant length is rho (standard error 0.005)
  expect_equal(abs(mean(exp(1i * x))), 0.4, tolerance = 0.05)
  expect_length(rcard(4, mu = 1:4, rho = 0.3), 4)
})

test_that("card recycles and rejects rho outside [0, 0.5]", {
  expect_equal(dcard(0.5, 0, c(0.1, 0.4)), c(dcard(0.5, 0, 0.1), dcard(0.5, 0, 0.4)))
  expect_error(dcard(1, 0, 0.6), "rho")
  expect_error(pcard(1, 0, -0.1), "rho")
  expect_error(rcard(1, 0, 0.7), "rho")
})
